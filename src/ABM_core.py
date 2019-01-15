#!/usr/bin/env python
# encoding: utf-8

"""
*****************************************************************************
Author(s):	Maximilian King  (previous authors: Lars Seemann - lseemann@uh.edu)
Email: Maximilian_King@brown.edu
Organization: Marshall Lab, Department of Epidemiology - Brown University

Description:
    Module responsible for ABM simulation events. Operates main loop over simulation run.
    Handles agent pairing, interaction, disease propagation, interventions, deaths, etc.


Copyright (c) 2016, Maximilian King
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************
"""

# Imports
import random
from random import Random
#from copy import deepcopy, copy
import os
import time
import collections


from scipy.stats import binom
from scipy.stats import poisson
from functools import wraps

try:
    from HIVABM_Population import PopulationClass, print_population
except ImportError:
    raise ImportError("Can't import PopulationClass")

try:
    from network_graph_tools import *
except ImportError, e:
    raise ImportError("Can't import network_graph_tools! %s" % str(e))

try:
    from ABM_partnering import *
except ImportError, e:
    raise ImportError("Can't import ABM_partnering! %s" % str(e))

try:
    from analysis_output import * #assessment_lib import *   #OLD FILE
except ImportError, e:
    raise ImportError("Can't import analysis_output! %s" % str(e))



PROF_DATA = {}
def profile(function):
    @wraps(function)
    def with_profiling(*args, **kwargs):
        start_time = time.time()

        ret = function(*args, **kwargs)

        elapsed_time = time.time() - start_time

        if function.__name__ not in PROF_DATA:
            PROF_DATA[function.__name__] = [0, []]
        PROF_DATA[function.__name__][0] += 1
        PROF_DATA[function.__name__][1].append(elapsed_time)
        #print elapsed_time

        return ret

    return with_profiling

def print_prof_data():
    for fname, data in PROF_DATA.items():
        max_time = max(data[1])
        avg_time = sum(data[1]) / len(data[1])
        print "Function %s called %d times. " % (fname, data[0])
        print '\tExecution time max: %.3f, average: %.3f, total %.3f' % (max_time, avg_time, sum(data[1]))

def clear_prof_data():
    global PROF_DATA
    PROF_DATA = {}


class HIVModel(NetworkClass):
    """
    :Purpose:
        This is the core class used to simulate
        the spread of HIV and drug use in one MSA
        (Metropolitan Statistical Area).

    :Input:
        N : int
            Number of agents. Default: 1000
        tmax: int
            Number of simulation steps (years).

        :py:class:`SocialNetworkClass` : Inherited
        :py:class:`PopulationClass` : Inherited

    :Attributes:
        :py:attr:`tmax` : int
            Number of time steps simulated.

        :py:attr:`CleanSyringeUsers` : list

        :py:attr:`SEPAgents` : dict
            Dictionary of users who participated in a
            syringe exchange program (SEP) {agent:time}.

        :py:attr:`DrugTreatmentAgents` : dict
            Dictionary of users who underwent drug
            treatment {agent:time}.

        :py:attr:`TestedAgents` : list
            List of agents who get tested for HIV every time step.

        :py:attr:`tmp_Agents` : dict
            Changes resulting from parsing through the agents
            and applying the update rules are stored
            in :py:attr:`tmp_agent_dict`.

        All attributes from :py:class:`SocialNetworkClass` \n

        All attributes from :py:class:`PopulationClass`

    :Methods:
        :py:meth:`_update_population` \n
        :py:meth:`_needle_transmission` \n
        :py:meth:`_sex_transmission` \n
        :py:meth:`_drug_transition` \n
        :py:meth:`_update_IDU` \n
        :py:meth:`_update_NIDU_ND` \n
        :py:meth:`_update_AllAgents` \n
        :py:meth:`_VCT` \n
        :py:meth:`_SEP` \n
        :py:meth:`_enter_drug_treatment` \n
        :py:meth:`_initiate_HAART` \n
        :py:meth:'_discontinue_HAART' \n
        :py:meth:`_get_partner` \n
        :py:meth:`_update_AllAgents` \n
        :py:meth:`run` \n
        :py:meth:`store_results` \n
        :py:meth:`get_HIV_prevalence_drugs` \n
        :py:meth:`get_HIV_prevalence` \n
        :py:meth:`plot_results` \n
        All methods from :py:class:`SocialNetworkClass` \n
        All methods from :py:class:`PopulationClass`
    """
    def __repr__(self):
        returnStr = "\n"
        returnStr += "Seed: %d\n"%(self.runseed)
        returnStr += "Npop: %d\n"%(params.N_POP)
        returnStr += "Time: %d\n"%(params.TIME_RANGE)
        returnStr += "Mode: %s\n"%(params.model)

        return returnStr

    def __init__(self, N, tmax, parameter_dict, runseed, popseed, netseed,
                 runtime_diffseed=False, model=None, network_type=None, HIVABM_Agent_set=None):
        """ Initialize HIVModel object """
        # Ensure param variable is are defined. For backwards compatibility with params.py files
        try:
            params.drawEdgeList
        except AttributeError:
            params.drawEdgeList = False

        try:
            params.inc_treat_HRsex_beh
        except AttributeError:
            params.inc_treat_HRsex_beh = False

        try:
            params.inc_treat_IDU_beh
        except AttributeError:
            params.inc_treat_IDU_beh = False

        try:
            params.calcNetworkStats
        except AttributeError:
            params.calcNetworkStats = False

        if (type(tmax) is not int):
            raise ValueError("Number of time steps must be integer")
        else:
            self.tmax = tmax

        if (type(runseed) is not int):
            raise ValueError("Random seed must be integer")
        elif runseed == 0:
            self.runseed = random.randint(1,1000000)
        else:
            self.runseed = runseed

        if (type(popseed) is not int):
            raise ValueError("Random seed must be integer")
        elif popseed == 0:
            self.popseed = random.randint(1,1000000)
        else:
            self.popseed = popseed

        if (type(netseed) is not int):
            raise ValueError("Random seed must be integer")
        elif popseed == 0:
            self.netseed = random.randint(1,1000000)
        else:
            self.netseed = netseed

        self.uniqueSeedID = 'r'+ str(runseed) + '_p' + str(popseed) + '_n' + str(netseed)



        self.current_dir = os.getcwd()
        print("=== Begin Initialization Protocol ===\n")
        self.ExistingLinksCollapsedList = list()

        # Computation runtime tic/tocs
        self.TIC = 0
        self.TOC = 0


        print("\tDictionary Read")

        # Risk network replaced social network
        if HIVABM_Agent_set:
            print("\tReading prefab agent set for population")
            self.All_agentSet = network_type
            self.Relationships = Agent_set(1,"Relationships")
        else:
            print("\tCreating population Class")
            #PopulationClass = network_type
            #PopulationClass.__init__(self, n=N, rSeed = rseed, model=model)
            
            #asdf.All_agentSet.remove_agent(asdf.All_agentSet.random_agent())
            #asdf.All_agentSet.print_subsets()

            #thisPopClass = PopulationClass(N, rseed, model)
            NetworkClass.__init__(self, N=N, network_type=network_type, popSeed=self.popseed, netSeed=self.netseed)
            self.All_agentSet.print_subsets()
            # thing.All_agentSet.print_subsets()
            # thing = PopulationClass(N, rseed, model)
            # thing.All_agentSet.print_subsets()
            # print type(thing.Trt_ART_agentSet)

        self.AdjMat = 0
        self.AdjMats_by_time = 0

        # keep track of current time step globally for dynnetwork report
        self.TimeStep = 0
        self.totalIncarcerated = 0


        print("\n\tCreating lists")
        # Other lists / dictionaries

        self.NewInfections = Agent_set(3, "NewInfections")
        self.NewDiagnosis = Agent_set(3, "NewDiagnosis")
        self.NewIncarRelease = Agent_set(3, "NewIncarRelease")
        self.NewHRrolls = Agent_set(3, "NewHRrolls")

        # Assess the distribution of number of interactions per timestep for each agent type
        self.ND_NumPartners = {'ND': [], 'NIDU': [], 'IDU': [], 'MSM': []}  # final counts
        self.NIDU_NumPartners = {'ND': [], 'NIDU': [], 'IDU': [], 'MSM': []}  # final counts
        self.IDU_NumPartners = {'ND': [], 'NIDU': [], 'IDU': [], 'MSM': []}  # final counts
        self.MSM_NumPartners = {'ND': [], 'NIDU': [], 'IDU': [], 'MSM': []}  # final counts

        self.tmp_ND_NumPartners_Count = {}
        self.tmp_NIDU_NumPartners_Count = {}
        self.tmp_IDU_NumPartners_Count = {}
        self.tmp_MSM_NumPartners_Count = {}
        self.tmp_WSW_NumPartners_Count = {}

        self.Acute_agents = []
        self.Transmit_from_agents = []
        self.Transmit_to_agents = []
        self.Transmission_tracker = {'SEX_MSM': {1: 0}, 'SEX_NMSM': {1: 0}, 'NEEDLE': {1: 0}}
        self.totalDiagnosis = 0
        self.treatmentEnrolled = False

        self.ResultDict = initiate_ResultDict()

        # Set seed format. 0: pure random, -1: Stepwise from 1 to nRuns, else: fixed value
        

        print "\tRun seed was set to:", runseed
        self.runRandom = Random(runseed)
        random.seed(self.runseed)
        np.random.seed(self.runseed)
        print "\tFIRST RANDOM CALL %d" %random.randint(0,100)

        print("\tReseting death count")
        self._reset_death_count()  # Number of death

        print("\tCreating network graph")
        self.create_graph_from_agents(self.All_agentSet)
        #self.get_Graph = NetworkClass(1000,m_0=1)

        print("\n === Initialization Protocol Finished ===")

    def run(self, save_adjlist_flag=1, dir_prefix='Results'):
        """
        Core of the model:
            1. Prints networkReport for first agents.
            2. Makes agents become HIV (used for current key_time tracking for acute)
            3. Loops over all time steps
                a. _update AllAgents()
                b. _reset_death_counts()
                c. _ self._die_and_replace()
                d. self._update_population()
                e. self._reset_partner_count()
        """
        def getStats(t):
            self.filler = 0
            print_stats(self, self.runseed, t
                ,self.All_agentSet
                ,self.HIV_agentSet
                ,self.incarcerated_agentSet
                ,self.Trt_PrEP_agentSet
                ,self.NewInfections
                ,self.NewDiagnosis
                ,self.num_Deaths
                ,self.ResultDict
                ,self.Relationships
                ,self.NewHRrolls
                ,self.NewIncarRelease)

        def print_components(t):
            name = 'componentReport_ALL'
            compReport = open('results/'+name+'.txt', 'a')
            components = sorted(nx.connected_component_subgraphs(self.G), key=len, reverse=True)
            compID = 0
            for comp in components:
                totN = nhiv = ntrtmt = ntrthiv = nprep = 0
                for ag in comp.nodes():
                    totN += 1
                    if (ag._HIV_bool):
                        nhiv += 1
                        if ag._treatment_bool:
                            ntrthiv += 1
                    elif (ag._treatment_bool):
                        ntrtmt += 1
                        if ag._PrEP_bool:
                            nprep += 1


                compReport.write("{rseed}\t{pseed}\t{nseed}\t{t}\t{compID}\t{totalN}\t{Nhiv}\t{Ntrtmt}\t{Nprep}\t{NtrtHIV}\n".format(
                    rseed=self.runseed,
                    pseed=self.popseed,
                    nseed=self.netseed,
                    t=t,
                    compID=compID,
                    totalN=totN,
                    Nhiv=nhiv,
                    Ntrtmt=ntrtmt,
                    Nprep=nprep,
                    NtrtHIV=ntrthiv))

                compID += 1
            compReport.close()

        #print "RANDOM CALL %d" %random.randint(0,100)


        def burnSimulation(burnDuration):
            print("\n === Burn Initiated for {} timesteps ===".format(burnDuration+1))
            for t in range(0, burnDuration + 1):
                # print '\n\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t.: BURN', t
                self._update_AllAgents(t, burn=True)

                if params.flag_DandR:
                    #print("\t\tdie and replace")
                    self._die_and_replace()

            # self.All_agentSet.print_subsets()
            print "\tBurn Cuml Inc:\t{}".format(self.NewInfections.num_members())
            self.NewInfections.clear_set()
            self.NewDiagnosis.clear_set()
            self.NewHRrolls.clear_set()
            self.NewIncarRelease.clear_set()
            # getStats(0)
            print(" === Simulation Burn Complete ===")

        # self.get_Graph = NetworkClass(1000,m_0=1)
        burnSimulation(params.burnDuration)

        # print "RANDOM CALL %d" %random.randint(0,100)


        print("\n === Begin Simulation Run ===")
        #print("\t Writing Agents to dynNet Report")
        if params.drawFigures:
            #self.get_Graph.draw_histogram(0)
            nNodes = self.G.number_of_nodes()
            self.visualize_network(coloring=params.drawFigureColor, 
                node_size=5000./nNodes, 
                curtime=0, 
                iterations=10, 
                label="Seed"+str(self.runseed))
        if params.calcComponentStats:
            print_components(0)
        # write agents to dynnetworkReport
        #self._writeDNR()

        # print "TOTAL AGENTS:", self.All_agentSet.num_members()
        self.cumInfT = 0
        self.cumInfW = 0
        self.cumInfB = 0

        def makeAgentZero(numPartners):
            firstHIV = self.runRandom.choice(self.DU_IDU_agentSet._members)
            i=0
            while i <= numPartners:
                    update_partner_assignments(self, 10000.0, self.get_Graph(), agent=firstHIV)
                    i += 1
            self._become_HIV(firstHIV, 0)
        #degree_sequence = sorted([d for n, d in self.get_Graph.G.degree()], reverse=True)
        #print degree_sequence
        #print firstHIV

        #self._become_HIV(firstHIV, 0)

        print("\t===! Start Main Loop !===")

        #If we are using an agent zero method, create agent zero.
        if params.flag_agentZero:
            makeAgentZero(4)
        #print(self.runtime_diffseed)
        #print(self.All_agentSet._subset)
        # if self.runtime_diffseed:
        #     print("setting rseed for post agent making")
        #     self.runRandom.seed()
        for t in range(1, self.tmax + 1):
            print '\n\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t.: TIME', t
            #print "RANDOM CALL %d" %random.randint(0,100)
            if params.drawFigures and t%params.intermPrintFreq == 0:
                #self.get_Graph.draw_histogram(0)
                self.visualize_network(coloring=params.drawFigureColor, 
                    node_size=5000./nNodes, 
                    curtime=t, 
                    iterations=10, 
                    label="Seed"+str(self.runseed))
                # self.visualize_network(coloring=params.drawFigureColor, 
                #     node_size=10, 
                #     curtime=t, 
                #     txtboxLabel=4, 
                #     iterations=10, 
                #     label=params.label)
            #todo: GET THIS TO THE NEW HIV COUNT
            print "\tSTARTING HIV count:%d\tTotal Incarcerated:%d\tHR+:%d\tPrEP:%d" % (self.HIV_agentSet.num_members(), self.incarcerated_agentSet.num_members(), self.highrisk_agentsSet.num_members(), self.Trt_PrEP_agentSet.num_members())
            #self.All_agentSet.print_agents()
            self.TimeStep = t

            #print self.All_agentSet._members
            self._update_AllAgents(t)

            #print "Results Dictionary update"
            getStats(t)

            #print "Reseting death count"
            self._reset_death_count()

            if params.flag_DandR:
                #print("\t\tdie and replace")
                self._die_and_replace()

            #print "\t\tENDING HIV count:%2.2f\tIncarcerated:%d\tHR+:%d"%(self.All_agentSet._subset["HIV"].num_members()/self.All_agentSet.num_members(), self.IncarceratedClass.num_members(),self.PrEP_agents_class.num_members()) #,self.HighriskClass.num_members())
            print "Number of relationships: %d"%self.Relationships.num_members()
            tested = len([tmpA for tmpA in self.HIV_agentSet._members if tmpA._tested])
            # print "Number tested: %d\t%.2f"%(tested, 1.0*tested/max(1,self.HIV_agentSet.num_members()))
            self.All_agentSet.print_subsets()

            newInfB = len([tmpA for tmpA in self.NewInfections._members if tmpA._race == 'BLACK'])
            newInfW = len([tmpA for tmpA in self.NewInfections._members if tmpA._race == 'WHITE'])
            newInfT = len(self.NewInfections._members)
            self.cumInfB += newInfB
            self.cumInfW += newInfW
            self.cumInfT += newInfT

            # print "\n\tGroup\tMo\tCuml"
            # print "\tTotal:\t%d\t%d"%(newInfT,self.cumInfT)
            # print "\tWhite:\t%d\t%d"%(newInfW,self.cumInfW)
            # print "\tBlack:\t%d\t%d"%(newInfB,self.cumInfB)
            # print self.Relationships.print_agent_relationshps()
            self.totalDiagnosis += len(self.NewDiagnosis._members)
            if self.totalDiagnosis > params.initTreatment and not self.treatmentEnrolled:
                self._enroll_treatment(t)


            self.NewInfections.clear_set()
            self.NewDiagnosis.clear_set()
            self.NewHRrolls.clear_set()
            self.NewIncarRelease.clear_set()
            self.num_Deaths


            #If set to draw the edge list, print list at each timestep
            if params.drawEdgeList:
                print "Drawing network edge list to file"
                fh=open("results/network/Edgelist_t{}.txt".format(t),'wb')
                self.write_G_edgelist(fh)
                fh.close()

            #self.get_Graph.draw_histogram()
            #print self.get_Graph.stat_connectivity()


        #self.get_Graph.visualize_network(iterations=5)
        #self.get_Graph.vizualize_network_graphviz(program='neato', coloring='Tested', time=t)
        #self.get_Graph.vizualize_network_graphviz(program='neato', coloring='SO', time=t)


        print_prof_data()
        #print params.PrEP_type
        #print params.PrEP_Target

        #print(self.All_agentSet._subset)
        if t%params.intermPrintFreq == 0:
            if params.calcNetworkStats:
                self.write_network_stats(t=t)
            if params.calcComponentStats:
                print_components(t)


    #@profile
    def _update_AllAgents(self, time, burn=False):
        """
        :Purpose:
            Update IDU agents:
            For each agent:
                1 - determine agent type
                2 - get partners
                3 - agent interacts with partners
                4 - drug transition
                5 - VCT (Voluntsry Counseling and Testing)
                6 - if IDU: SEP, treatment
                7 - if HIV: HAART, AIDS
                8 - drug cessation

        :Input:
            agent, time

        :Output:
            none
®
        """
        num_HIV = len(self.HIV_agents)
        #print("\t\t= Begin Agents Partnering =")
        if time == 0:
            i=0
            # self.create_graph_from_agents(self.All_agentSet)

            # self.get_Graph.plot_DegreeDistribution()
            #self.vizualize_network_graphviz(program='neato', coloring='Tested')
            # self.visualize_network(node_size=10, coloring="HIV")
            # self.visualize_network(node_size=10, coloring="Trtmt")
            # update_partner_assignments(self, params.PARTNERTURNOVER, self.get_Graph)
            #self.get_Graph.create_graph_from_relationships(self.Relationships)
        elif params.flag_staticN == False:
            update_partner_assignments(self, params.PARTNERTURNOVER, self.get_Graph)

        else:
            pass
            #self.get_Graph.vizualize_network_graphviz(program='neato', coloring='Tested', time=time)
            # if time%12==0:self.get_Graph.plot_DegreeDistribution(time)
        
        #print("\t\t= Updated Partners =")
        self.Acute_agents = []
        self.Transmit_from_agents = []
        self.Transmit_to_agents = []
        
        #print("\t\t= Printed DynNetwork Report =")
        #print("\t\t= Begin Agents Operations =")
        #print("\t\t= Relationship Iterations =")
        #print("\n\nSTARTING RELATIONSHIPS")
        # self.Relationships.print_agent_relationshps()
        
        for rel in self.Relationships._members:
            #print "Rel: ",rel
            if burn:
                pass
            else:
                self._agents_interact(rel._ID1, rel._ID2, time, rel)
            if params.flag_staticN:
                pass
            else:
                if rel.progress():
                    try:self.get_Graph().remove_edge(rel._ID1, rel._ID2)
                    except:pass
                    self.Relationships.remove_agent(rel)
                    # relID = self.Relationships._members.index(rel)
                    # self.Relationships._members.pop(relID)
                    #print self.Relationships.is_member(rel)

                    #print rel
                    #self.Relationships._members.discard(rel)
                    del rel

                #pass
        # print("\n\nENDING RELATIONSHIPS")
        # print type(self.Relationships)
        # print self.Relationships._members
            #print self.Relationships.num_members()

        if params.flag_HR:
            #print("\t\t= High Risk Group functions =")
            for tmpA in self.highrisk_agentsSet.iter_agents():
                if tmpA._highrisk_time > 0:
                    tmpA._highrisk_time -= 1
                else:
                    self.highrisk_agentsSet.remove_agent(tmpA)
                    tmpA._highrisk_bool = False
                    
                    if params.model == 'Incar':
                        if tmpA._SO == "HM":
                            tmpA._mean_num_partners -= params.HR_partnerScale
                        elif tmpA._SO == "HF":
                            tmpA._mean_num_partners -= params.HR_partnerScale

        #print("\t\t= Agents Iterations (Incar/test/AIDS/HAART/PrEP =")
        
        #random.shuffle(self.All_agentSet._members)
        for agent in self.All_agentSet.iter_agents():#self.Agents: #NEW METHOD
            
            agent_drug_type = agent._DU
            agent_sex_type = agent._SO
            agent_HIV_status = agent._HIV_bool

            agent_incarcerated = agent._incar_bool#agent_dict['incar_t']

            #print("\tEnded partner interactions\n")
            #self._VCT(agent, time)  # VCT: Counseling and Testing
            #print "\tIncarcerate Agents"
            agent._timeAlive += 1
            if params.flag_incar:# and not burn:
                #print("\n\t\t= Incarceration effects =")
                self._incarcerate(agent, time)

            #if agent_drug_type == 'IDU':
                #self._SEP(agent, time)  # SEP: Syringe Exchange program
            
            if agent_drug_type in ['NIDU', 'IDU'] and False:
                #print("\tDrug Cessation")
                self._drug_cessation(agent, agent_drug_type)
                #print("\tEnter/Exit Drug Treatment")
                self._enter_and_exit_drug_treatment(agent, time)
            

            if agent_HIV_status:
                if burn:
                    #print "Viral Load Reset"
                    #self._become_HIV(agent, time)  # reset viral load ############# TURNED OFF NO NEED TO RESET VIRAL LOAD
                    #P_HAART = self.runRandom.random()
                    #if P_HAART < 0.005:
                    #print("\n\t\t= Testing agents =")
                    if agent._incar_treatment_time >= 1:
                        agent._incar_treatment_time -= 1


                self._HIVtest(agent, time)
                self._progress_to_AIDS(agent, agent_drug_type)

                if params.flag_ART:
                    self._initiate_HAART(agent, time)
                    #hiv_t = agent_dict['HIV_t']
                    #hiv_t += 1
                    agent._HIV_time += 1
                    # self.tmp_Agents[agent].update({'HIV_t': hiv_t})
            else:
                if params.flag_PrEP:
                    if time >= params.PrEP_startT:
                        if agent._PrEP_bool:
                            self._discont_PrEP(agent, time)
                        elif params.PrEP_target_model == 'Clinical':
                            pass
                        elif params.PrEP_target_model == 'RandomTrial':
                            pass
                        elif self._PrEP_elligible(agent, time) and not agent._PrEP_bool:
                            self._initiate_PrEP(agent, time)
                        # if not agent._PrEP_bool:
                        #     self._initiate_PrEP(agent, time)
        #print("\t\t= End Agents Operations =")


        if params.flag_PrEP and time >= params.PrEP_startT:
            if params.PrEP_target_model == 'Clinical':
                if time > params.PrEP_startT:
                    numPrEP_agents = self.Trt_PrEP_agentSet.num_members()
                    target_PrEP = int((self.All_agentSet.num_members()-self.All_agentSet._subset["HIV"].num_members()) * params.PrEP_Target)
                    elligiblePool = [ag for ag in self.All_agentSet._subset['SO']._subset['MSM']._members if (ag._PrEP_bool == False and ag._HIV_bool == False)]
                    # print "Eligible PrEP pool size: ",len(elligiblePool)

                    while(numPrEP_agents < target_PrEP):
                        numPrEP_agents = self.Trt_PrEP_agentSet.num_members()
                        target_PrEP = int((self.All_agentSet.num_members()-self.All_agentSet._subset["HIV"].num_members()) * params.PrEP_Target)
                        # print "%d/%d"%(numPrEP_agents, target_PrEP)
                        self._initiate_PrEP(self._get_clinic_agent(params.PrEP_clinic_cat, elligiblePool), time)
                # elif self._PrEP_elligible(agent, time):
                #     self._initiate_PrEP(agent, time)
            elif params.PrEP_target_model == 'RandomTrial' and time == params.PrEP_startT:
                print "Starting random trial"
                components = sorted(nx.connected_component_subgraphs(self.G), key=len, reverse=True)
                totNods = 0
                for comp in components:
                    totNods += comp.number_of_nodes()
                    if self.runRandom.random() < 0.5:
                        #Component selected as treatment pod!
                        for ag in comp.nodes():
                            if (ag._HIV_bool == False) and (ag._PrEP_bool == False):
                                ag._treatment_bool = True
                                if self.runRandom.random() < params.PrEP_Target:
                                    self._initiate_PrEP(ag, time, force=True)
                print "Total agents in trial: ",totNods
        #print("\t\t !!!! ALL AGENTS UPDATED !!!\n")


    def _agents_interact(self, agent_1, agent_2, time, rel):
        """
        :Purpose:
            Let IDU agent interact with a partner.
            Update IDU agents:
                1 - determine transition type
                2 - Injection rules
                3 - Sex rules
                4 - HIV transmission
                5 - SEP

        :Input:
            agent : int

            partner : int

            time : int

        Output:
            none

        """
        #print agent
        partner_HIV_status = agent_2._HIV_bool#self.get_agent_characteristic(partner, 'HIV')
        agent_HIV_status = agent_1._HIV_bool#self.get_agent_characteristic(agent, 'HIV')
        agent_incar = agent_1._incar_bool
        partner_incar = agent_2._incar_bool
        eligible = False


        #If either agent is incarcerated, skip their interaction
        if agent_incar or partner_incar:
            return

        #Else if neither agent is HIV (shouldn't be possible), skip their interaction to save computation time
        elif not agent_HIV_status and not partner_HIV_status:
            #print "\t\t!!! Neither agents HIV+, moving on A:%d P:%d"%(agent_1.get_ID(), agent_2.get_ID())
            return

        elif agent_HIV_status: #If agent_1 is HIV
            if partner_HIV_status: #If agent_1 and agent_2 are both HIV, skip interaction
                #print "\t\t!!! BOTH agents HIV+, moving on A:%d P:%d"%(agent_1.get_ID(), agent_2.get_ID())
                return
            else:#Agent is HIV, partner is succept
                agent = agent_1
                partner = agent_2
                eligible = True
        elif partner_HIV_status: #If agent_2 is HIV and we have tested both HIV +/-, agent_2 is HIV, agent_1 is succept
            agent = agent_2
            partner = agent_1
            eligible = True


        if eligible:
            #print("Elligible\t%d\t%d"%(agent._ID, partner._ID))
            partner_drug_type = partner._DU  # self.get_agent_characteristic(partner, 'Drug Type')
            agent_drug_type = agent._DU  # self.get_agent_characteristic(agent, 'Drug Type')
            partner_sex_type = partner._SO  # self.get_agent_characteristic(partner, 'Sex Type')
            agent_sex_type = agent._SO  # self.get_agent_characteristic(agent, 'Sex Type')
            partner_HIV_status = partner._HIV_bool  # self.get_agent_characteristic(partner, 'HIV')
            agent_HIV_status = agent._HIV_bool  # self.get_agent_characteristic(agent, 'HIV')
            agent_incar = agent._incar_bool
            partner_incar = partner._incar_bool
            if partner_drug_type == 'IDU' and agent_drug_type == 'IDU':
                # Injection is possible
                #If agent is on post incar HR treatment to prevent IDU behavior, pass IUD infections
                if agent._incar_treatment_time > 0 and params.inc_treat_IDU_beh:
                    pass

                elif self._sex_possible(agent_sex_type, partner_sex_type):
                    # Sex is possible
                    rv = self.runRandom.random()
                    if rv < 0.6: #Needle only (60%)
                        #print "Needle inc (IDUs)"
                        self._needle_transmission(agent, partner, time)
                    elif rv < 0.6 + 0.2: #Sex only (20%)
                        #print "Sex inc (IDUs)"
                        self._sex_transmission(agent, partner, time, rel)  # , num_interactions)
                    else: #Both sex and needle (20%)
                        #print "Needle and sex inc (IDUs)"
                        self._needle_transmission(agent, partner, time)
                        self._sex_transmission(agent, partner, time, rel)  # , num_interactions)
                else:
                    # Sex not possible, needle only
                    #print "Needle inc (IDUs)"
                    self._needle_transmission(agent, partner, time)

            elif (partner_drug_type in ['NIDU', 'NDU'] or agent_drug_type in ['NIDU', 'NDU']):
                #print "Sex inc (ND/NIDU)"
                if self._sex_possible(agent_sex_type, partner_sex_type):
                    self._sex_transmission(agent, partner, time, rel)  # ,num_interactions)
                else:
                    return
                    #print "!!!!!!!!!!!!SEX NOT POSSIBLE A:%s \tP:%s"%(agent_sex_type, partner_sex_type)
            else:
                raise ValueError("Agents must be either IDU, NIDU, or ND")



    def _drug_transition(self, agent, partner):
        """
        :Purpose:
            Simulate transition of drug behavior. The following scenarios are
            possible:
            + ND agent might become NIDU when meeting NIDU
            + NIDU might become IDU when meeting IDU
            The function is only applied for NIDU and ND users.

        :Input:
            agents : int
            partner : int

        :Output: -
        """

        partner_drug_type = self.get_agent_characteristic(partner, 'Drug Type')
        agent_drug_type = self.get_agent_characteristic(agent, 'Drug Type')
        Flag_Partner_IDU_NIDU_Transition = 0
        Flag_Agent_IDU_NIDU_Transition = 0

        # NIDU -> IDU
        if agent_drug_type == 'NIDU' and partner_drug_type == 'IDU':
            if self.runRandom.random() < 0.00875 / 12:
                #print "Agent %d just became IDU" % (agent, )
                self.tmp_Agents[agent].update({'Drug Type': 'IDU'})  # agent becomes IDU
            # Sex type lists
            if agent in self.tmp_NIDU_agents:
                self.tmp_NIDU_agents.remove(agent)
            if agent in self.tmp_ND_agents:  # agent might have transitioned into ND before
                self.tmp_ND_agents.remove(agent)
            if agent not in self.tmp_IDU_agents:
                self.tmp_IDU_agents.append(agent)

        elif partner_drug_type == 'NIDU' and agent_drug_type == 'IDU':
            if self.runRandom.random() < 0.0175 / 12:
                self.tmp_Agents[partner].update({'Drug Type': 'IDU'})  # partner becomes IDU
            # Sex type lists
            if partner in self.tmp_NIDU_agents:
                self.tmp_NIDU_agents.remove(partner)
            if partner in self.tmp_ND_agents:  # agent might have transitioned into ND before
                self.tmp_ND_agents.remove(partner)
            if partner not in self.tmp_IDU_agents:
                self.tmp_IDU_agents.append(partner)

        ## ND -> IDU
        elif agent_drug_type == 'ND' and partner_drug_type == 'IDU':
            if self.runRandom.random() < 0.001:
                self.tmp_Agents[agent].update({'Drug Type': 'IDU'})  # agent becomes NIDU
            # Sex type lists
            if agent in self.tmp_ND_agents:  # agent might have transitioned already
                self.tmp_ND_agents.remove(agent)
            if agent not in self.tmp_IDU_agents:
                self.tmp_IDU_agents.append(agent)

        # ND -> IDU
        elif partner_drug_type == 'ND' and agent_drug_type == 'IDU':
            if self.runRandom.random() < 0.001:
                self.tmp_Agents[agent].update({'Drug Type': 'IDU'})  # agent becomes NIDU
            # Sex type lists
            if agent in self.tmp_ND_agents:  # agent might have transitioned already
                self.tmp_ND_agents.remove(agent)
            if agent not in self.tmp_IDU_agents:
                self.tmp_IDU_agents.append(agent)

        # ND -> NIDU
        elif agent_drug_type == 'ND' and partner_drug_type == 'NIDU':
            if self.runRandom.random() < 0.005:
                self.tmp_Agents[agent].update({'Drug Type': 'NIDU'})  # agent becomes NIDU
            # Sex type lists
            if agent in self.tmp_ND_agents:  # agent might have transitioned already
                self.tmp_ND_agents.remove(agent)
            if agent not in self.tmp_NIDU_agents:
                self.tmp_NIDU_agents.append(agent)
        # ND -> NIDU
        elif partner_drug_type == 'ND' and agent_drug_type == 'NIDU':
            if self.runRandom.random() < 0.005:
                self.tmp_Agents[partner].update({'Drug Type': 'NIDU'})  # partner becomes NIDU
            # Sex type lists
            if partner in self.tmp_ND_agents:  # agent might have transitioned already
                self.tmp_ND_agents.remove(partner)
            if partner not in self.tmp_NIDU_agents:
                self.tmp_NIDU_agents.append(partner)
                # if partner in self.tmp_IDU_agents:     # agent might have previously transitioned into IDU
            #    self.tmp_IDU_agents.remove(agent)

        # NIDU -> ND from agent's perspective
        elif agent_drug_type == 'NIDU' and partner_drug_type == 'ND':
            if self.runRandom.random() < 0.001:
                self.tmp_Agents[agent].update({'Drug Type': 'ND'})
            # Sex type lists
            if agent in self.tmp_NIDU_agents:  # agent might have transitioned already
                self.tmp_NIDU_agents.remove(agent)
            if agent in self.tmp_IDU_agents:  # agent might have previously transitioned into IDU
                self.tmp_IDU_agents.remove(agent)
            if agent not in self.tmp_ND_agents:
                self.tmp_ND_agents.append(agent)
            if agent in self.DrugTreatmentAgents_current:
                self._exit_drug_treatment(agent)

        # NIDU -> ND from partner's perspective
        elif partner_drug_type == 'NIDU' and agent_drug_type == 'ND':
            if self.runRandom.random() < 0.001:
                self.tmp_Agents[partner].update({'Drug Type': 'ND'})
            # Sex type lists
            if partner in self.tmp_NIDU_agents:  # partner might have transitioned already
                self.tmp_NIDU_agents.remove(partner)
            if partner in self.tmp_IDU_agents:  # partner might have previously transitioned into IDU
                self.tmp_IDU_agents.remove(partner)
            if partner not in self.tmp_ND_agents:
                self.tmp_ND_agents.append(partner)
            if partner in self.DrugTreatmentAgents_current:
                self._exit_drug_treatment(partner)
        else:
            pass  # transition not possible


    def get_acute_status(self, agent, time):
        """
        :Purpose:
            Simulate random transmission of HIV between two IDU agents
            through needle.\n
            Needed in _update_IDUand
        :Input:
            agents : int
            partner : int
        time : int
        :Output: -
        """
        acuteTimePeriod = 3
        hiv_t = agent._HIV_time#self.Agents[agent]['HIV_t']



        #if time > 0:
            #print self.HIV_key_transitiontime
        if hiv_t <= acuteTimePeriod and hiv_t > 0:
            #print "Agent %d has been sick for %d timesteps"%(agent, hiv_t)
            return True
        else:
            return False


    def get_transmission_probability(self, agent, interaction):
        """ Decriptor
            :Purpose:
            Determines the probability of a transmission event based on type. Determines if act is needle/sexual,

            :Input:
                N : int
                Number of agents. Default: 1000
                tmax: int
                Number of simulation steps (years).

                :py:class:`SocialNetworkClass` : Inherited
                :py:class:`PopulationClass` : Inherited

            :Attributes:
                :py:attr:`tmax` : int
                    Number of time steps simulated.
                """

        sex_type = agent._SO#self.get_agent_characteristic(agent, 'Sex Type')
        race_type = agent._race#self.get_agent_characteristic(agent, 'Race')
        #time = self.TimeStep
        ageBin = agent._ageBin
        tested = agent._tested#self.get_agent_characteristic(agent, 'Tested')
        onHAART = agent._HAART_bool
        # viral_load = self.Viral_load[agent]
        # v = 10 ** viral_load
        # p = (1. - (1. - (0.317 * (v ** 1.02)) / (v ** 1.02 + 13938 ** 1.02)) ** (1. / 83.17544)) * (0.014 / 0.003)

        agentAdherence = agent._HAART_adh#str(self.AdherenceAgents[agent])
        "Logic for if needle or sex type interaction"
        if interaction == 'NEEDLE':
            p = params.TransmissionProbabilities['NEEDLE'][str(agentAdherence)]

        elif interaction == 'SEX':
            p = params.TransmissionProbabilities['SEX'][sex_type][str(agentAdherence)]

        isAcute = self.get_acute_status(agent, 0)

        #Scaling parameter for acute HIV infections
        if(isAcute):
            p = p * params.cal_AcuteScaling

        #Scaling parameter for positively identified HIV agents
        if tested:
            p = p * (1 - params.cal_RR_Dx)

        #Tuning parameter for ART efficiency
        if onHAART:#self.AdherenceAgents[agent] > 0:
            #print "Agent on HAART, no xmission"
            p = p * params.cal_RR_HAART

        #Racial calibration parameter to attain proper race incidence disparity
        if race_type == 'BLACK':
            p = p * params.cal_raceXmission

        #Scaling parameter for per act transmission.
        p = p * params.cal_pXmissionScaling

        #Scaling parameter for age bin to attain proper age incidence disparity
        #p = p * params.cal_ageXmission[ageBin]


        return p

    def _needle_transmission(self, agent, partner, time):
        """
        :Purpose:
            Simulate random transmission of HIV between two IDU agents
            through needle.\n
            Needed in _update_IDUand
        :Input:
            agents : int
            partner : int
        time : int
        :Output: -
        """

        #Param to scale number of partners
        #NEEDLESCALINGPARAM = 0.7

        # both must be IDU
        partner_drug_type = partner._DU#self.get_agent_characteristic(partner, 'Drug Type')
        agent_drug_type = agent._DU#self.get_agent_characteristic(agent, 'Drug Type')
        agent_race = agent._race#self.get_agent_characteristic(agent, 'Race')
        agent_sex_type = agent._SO#self.get_agent_characteristic(agent, 'Sex Type')
        Race_Agent = agent._race#self.get_agent_characteristic(agent, 'Race')
        Type_agent = agent._SO#self.get_agent_characteristic(agent, 'Sex Type')


        if not (partner_drug_type == 'IDU' and agent_drug_type == 'IDU'):
            raise ValueError("To share a needle both agents must be IDU!%s %s" %
                             (str(agent_drug_type), str(partner_drug_type)))
        NumberP = len(self.ExistingLinksCollapsedList)
        # Do they share a needle?
        # OLD: if (agent in self.SEPAgents or partner in self.SEPAgents):
        #SEPstat = self._SEP(agent, time)
        SEPstat = agent._SNE_bool

        """
        tmpRand = self.runRandom.random()
        if (tmpRand > 0.43 and agent_race == 'WHITE') or (tmpRand > 0.27 and agent_race == 'BLACK'):################ FIX FOR RACIAL DISPARITY
            #print "SEP SAVE %s %d rolled %.3lf"%(agent_race, agent, tmpRand)
            SEPstat = True

        """

        isAcute = self.get_acute_status(agent, time)
        HIV_agent = agent._HIV_bool#self.get_agent_characteristic(agent, 'HIV')
        HIV_partner = partner._HIV_bool#self.get_agent_characteristic(partner, 'HIV')
        MEAN_N_ACTS = params.DemographicParams[Race_Agent][Type_agent]['NUMSexActs'] * params.cal_NeedleActScaling
        share_acts = poisson.rvs(MEAN_N_ACTS, size=1)

        if SEPstat:
            #print "SEP SAVE"
            p_UnsafeNeedleShare = 0.02 # no needle sharing
        else:  # they do share a needle
            # HIV+ ?
            
            #share_acts = int(random.uniform(1,30))
            if share_acts < 1:
                share_acts = 1

            p_UnsafeNeedleShare = params.DemographicParams[agent_race][agent_sex_type]['NEEDLESH'] * params.safeNeedleExchangePrev
            #MSexActs = self.ProbTables[Race_Agent][Type_agent]['NUMSexActs']

        for n in range(share_acts):
            if self.runRandom.random() > p_UnsafeNeedleShare:
                share_acts -= 1

        if HIV_agent == 1 and HIV_partner == 0 and share_acts >= 1.0:
            p = self.get_transmission_probability(agent, 'NEEDLE')
            #print p
            p_transmission = binom.pmf(1.0, share_acts, p)

            p_total_transmission = 0
            if share_acts == 1:
                p_total_transmission = p
            else:
                p_total_transmission = 1. - binom.pmf(0, share_acts, p)
                # for k in range(1, share_acts+1):
                #     temp = binom.pmf(k, share_acts, p)
                #     #print temp
                #     p_total_transmission += temp


            #print "\t\t\tHIV+ NED Act A:%d on P:%d\t Must be less than %.10lf\t(k:1 n:%.2lf, p:%.5lf) **OLD pTrans:%.5lf\tAcute:%s" % (agent, partner, p_total_transmission, share_acts, p, p_transmission, isAcute)
            if self.runRandom.random() < p_total_transmission:
                # if agent HIV+ partner becomes HIV+
                #transmit_HIV(agent, partners, t, )
                self._become_HIV(partner, time)
                #self.Transmission_tracker['NEEDLE'][time] += 1
                self.Transmit_from_agents += [agent]
                self.Transmit_to_agents += [partner]
                #if agent in ((self.HIV_key_transitiontime[time - 1] if time > 1 else [])
                #                 + (self.HIV_key_transitiontime[time - 2] if time > 2 else [])
                #                 + (self.HIV_key_transitiontime[time - 3] if time > 3 else [])):
                if(isAcute):
                    self.Acute_agents += [agent]
                    #print "\t\t(ACUTE)\tNE_HIV from agent %d to partner %d \t@ p=%.5lf transmissionp=%.5lf n:%d" %(agent._ID, partner._ID, p, p_total_transmission, share_acts)
                else:
                    pass
                    #print "\t\t\t\tNE_HIV from agent %d to partner %d \t@ p=%.5lf transmissionp=%.5lf n:%d" %(agent._ID, partner._ID, p, p_total_transmission, share_acts)

    #@profile
    def _sex_transmission(self, agent, partner, time, rel):

        """
        :Purpose:
            Simulate random transmission of HIV between two agents through Sex.
            Needed for all users. Sex is not possible in case the agent and
            assigned partner have incompatible Sex behavior.

        :Input:
            agents : int
            partner : int
            time : int
        number_of_interaction : int

        :Output:
            none
        """
        #SEXSCALINGPARAM = 0.2

        # Double check: Sex possible?
        Type_agent = agent._SO#self.get_agent_characteristic(agent, 'Sex Type')
        Type_partner = partner._SO#self.get_agent_characteristic(partner, 'Sex Type')
        if not self._sex_possible(Type_agent, Type_partner):
            raise ValueError("Sex must be possible! %s %s" % (
                str(Type_agent), str(Type_partner)))
        #NumberP = len(self.ExistingLinksCollapsedList)

        # HIV status of agent and partner
        # Everything from here is only run if one of them is HIV+

        HIVstatus_Agent = agent._HIV_bool#self.get_agent_characteristic(agent, 'HIV')
        HIVstatus_Partner = partner._HIV_bool#self.get_agent_characteristic(partner, 'HIV')
        AIDSstatus_Agent = agent._AIDS_bool#self.get_agent_characteristic(agent, 'AIDS')
        AIDSstatus_Partner = partner._AIDS_bool#self.get_agent_characteristic(partner, 'AIDS')
        Race_Agent = agent._race#self.get_agent_characteristic(agent, 'Race')
        isAcute = self.get_acute_status(agent, time)

        if HIVstatus_Partner: pass
        if HIVstatus_Agent and HIVstatus_Partner:
            #print "\t\t!!A%d %s P%d %s"%(agent.get_ID(), agent._HIV_bool, partner.get_ID(), partner._HIV_bool)
            #print "BOTH AGENTS %d and %d WERE HIV+!?!?!?" %(agent._ID, partner._ID)
            
            return
            #exit(10)
        elif HIVstatus_Agent == 1 or HIVstatus_Partner == 1:
            # Sex between men?
            if Type_agent == 'MSM' and Type_partner == 'MSM':
                SexBetweenMen = 1
            else:
                SexBetweenMen = 0
            # Define probabilities for unsafe sex

            # unprotected sex probabilities for primary partnerships
            p_UnsafeSafeSex1 = params.DemographicParams[Race_Agent][Type_agent]['UNSAFESEX']
            MSexActs = self._get_number_of_sexActs(agent) * params.cal_SexualActScaling
            #MSexActs = self.ProbTables[Race_Agent][Type_agent]['NUMSexActs'] * self.SEXSCALINGPARAM
            #print "Unsafe:%.5lf\tMSexActs:%.2lf\tOLDMSexActs:%.2lf"%(p_UnsafeSafeSex1,MSexActs,self.MEAN_S_ACTS)
            #print "MSEX",MSexActs
            T_sex_acts1 = int(poisson.rvs(MSexActs, size=1))

            num_int = rel._total_sex_acts
            #Get condom usage
            if num_int < 10:
                if num_int == 0:
                    p_UnsafeSafeSex1 = 0.443
                elif num_int == 1:
                    p_UnsafeSafeSex1 = 0.481
                else:
                    p_UnsafeSafeSex1 = 0.514
            else: #More than 10 acts
                p_UnsafeSafeSex1 = 0.759


            #Reduction of risk acts between partners for condom usage
            U_sex_acts1 = T_sex_acts1
            for n in range(U_sex_acts1):
                if self.runRandom.random() < p_UnsafeSafeSex1:
                    U_sex_acts1 -= 1


            U_sex_acts2 = U_sex_acts1
            #Reduction of risk acts between partners for PrEP adherence
            # for n in range(U_sex_acts1):
            #     if agent._PrEP_bool or partner._PrEP_bool:
            #         if agent._PrEP_adh == 1 or partner._PrEP_adh == 1:
            #             if self.runRandom.random() < 0.96:
            #                 U_sex_acts2 -= 1
            #         else:
            #             if self.runRandom.random() < 0.76:
            #                 U_sex_acts2 -= 1


            #U_sex_acts2 = 1
            #print "MeanS_act:%.2lf\tT_sex_acts1:%.2lf\tp_UnsafeSex1:%.2lf\tU_sex_acts1:%.2lf"%(MSexActs, T_sex_acts1,p_UnsafeSafeSex1,U_sex_acts1)
            if U_sex_acts2 >= 1:
                # if agent HIV+
                rel._total_sex_acts += U_sex_acts2
                if HIVstatus_Agent == 1 or HIVstatus_Partner == 1:
                    ppAct = self.get_transmission_probability(agent, 'SEX')

                    #Reduction of transmissibility for acts between partners for PrEP adherence
                    if agent._PrEP_bool or partner._PrEP_bool:
                        if agent._PrEPresistance or partner._PrEPresistance:
                            pass

                        elif params.PrEP_type == 'Oral':
                            if agent._PrEP_adh == 1 or partner._PrEP_adh == 1:
                                ppAct = ppAct * (1.0-params.PrEP_AdhEffic) #0.04
                            else:
                                ppAct = ppAct * (1.0-params.PrEP_NonAdhEffic) #0.24

                        elif params.PrEP_type == 'Inj':
                            ppActReduction = -1.0*np.exp(-5.528636721*partner._PrEP_load) + 1
                            #print "adwadw"
                            #print ppActReduction, partner._PrEP_load, partner._PrEP_lastDose
                            if agent._PrEP_adh == 1 or partner._PrEP_adh == 1:
                                ppAct = ppAct * (1.0-ppActReduction) #0.04



                    #p_transmission = binom.pmf(1, U_sex_acts1, p)

                    p_total_transmission = 0
                    if U_sex_acts2 == 1:
                        p_total_transmission = ppAct
                    else:
                        p_total_transmission = 1. - binom.pmf(0, U_sex_acts1, ppAct)
                        # for k in range(1, U_sex_acts2):
                        #     temp = binom.pmf(k, U_sex_acts1, ppAct)
                        #     #print temp
                        #     p_total_transmission += temp

                    #print "\t\t\tHIV+ SEX Act A:%d on P:%d\t Must be less than %.10lf\t(k:1 n:%.2lf, p:%.5lf) **OLD pTrans:%.5lf\tAcute:%s" % (agent, partner, p_total_transmission, U_sex_acts1, p, p_transmission, isAcute)
                    #print "\t\t\tHIV+ SEX ACT\tMust be less than %.10lf" % p_total_transmission
                    if self.runRandom.random() < p_total_transmission:


                        # if agent HIV+ partner becomes HIV+
                        self.Transmit_from_agents += [agent]
                        self.Transmit_to_agents += [partner]
                        #if Type_agent == 'MSM': self.Transmission_tracker['SEX_MSM'][time] += 1
                        #if Type_agent != 'MSM': self.Transmission_tracker['SEX_NMSM'][time] += 1

                        #print "\t\t\t\tST_HIV from agent %d to partner %d @ %.5lf" %(agent, partner, p)
                        self._become_HIV(partner, time)
                        #print 'INFECTION', Type_agent, DrugType_Agent
                        if(isAcute):
                            self.Acute_agents += [agent]
                            #print "\t\t(ACUTE)\tST_HIV from agent %d to partner %d \t@ p=%.5lf transmissionp=%.5lf n:%d" %(agent._ID, partner._ID, p, p_total_transmission, U_sex_acts1)
                        else:
                            pass
                        #print "\t\t\t\tST_HIV from agent %d to partner %d \t@ p=%.5lf transmissionp=%.5lf n:%d" %(agent._ID, partner._ID, ppAct, p_total_transmission, U_sex_acts1)


                        """if agent in ((self.HIV_key_transitiontime[time - 1] if time > 1 else [])
                                         or (self.HIV_key_transitiontime[time - 2] if time > 2 else [])
                                         or (self.HIV_key_transitiontime[time - 3] if time > 3 else [])):
                            self.Acute_agents += [agent]
                            print '\t\t\t\tACUTE', agent"""

                # if partner HIV+
                """if HIVstatus_Agent == 0 and HIVstatus_Partner == 1:
                    p = self.get_transmission_probability(partner, 'SEX')
                    p_transmission = binom.pmf(1, U_sex_acts1, p)
                    #print "\t\t\tHIV+ SEX ACT\tMust be less than %.10lf" % p_transmission
                    if self.runRandom.random() < p_transmission:
                        # if agent HIV+ partner becomes HIV+
                        self.Transmit_from_agents += [partner]
                        if Type_agent == 'MSM': self.Transmission_tracker['SEX_MSM'][time] += 1
                        if Type_agent != 'MSM': self.Transmission_tracker['SEX_NMSM'][time] += 1

                        print "\t\t\t\tST_HIV from partner %d to agent %d @ %.5lf" %(partner, agent, p)
                        self._become_HIV(agent, time)
                        if partner in ((self.HIV_key_transitiontime[time - 1] if time > 1 else [])
                                           or (self.HIV_key_transitiontime[time - 2] if time > 2 else [])
                                           or (self.HIV_key_transitiontime[time - 3] if time > 3 else [])):
                            self.Acute_agents += [partner]
                            print '\t\t\t\tACUTE', partner"""
            else:
                return
                #print "U_sex = %d < 1" % U_sex_acts1

    def _get_number_of_sexActs(self, agent):
        """
        :Purpose:
            agent becomes HIV agent. Update all appropriate list and
            dictionaries.

        :Input:
            agent : int

        """
        # 1 time per year 96 1.9 29 0.9 67 3.4
        # 2–5 times per year 428 8.2 184 5.8 244 12.2
        # 6–11 times per year 328 6.3 183 5.7 145 7.3
        # 12–23 times per year 376 7.2 251 7.9 125 6.3
        # 24–35 times per year 1,551 29.9 648 20.3 903 45.3
        # 36–51 times per year 1,037 20.0 668 20.9 369 18.5
        # 52–155 times per year 644 12.4 605 18.9 39 2.0
        # >156 times per year 733 14.1 631 19.7 102 5.1
        rv = self.runRandom.random()
        pMatch = 0.0
        i=0

        while(True):
            i += 1
            pMatch += params.sexualFrequency[i]['p_value']
            if rv <= pMatch:
                minSA = params.sexualFrequency[i]['min']
                maxSA = params.sexualFrequency[i]['max']
                return self.runRandom.randrange(minSA,maxSA,1)
            if i==5:break

        # if diceroll < 0.018:
        #     sexActs = 1
        # elif diceroll < 0.019 + 0.082:
        #     sexActs = self.runRandom.randrange(2,5,1)
        # elif diceroll < 0.019 + 0.082 + 0.063:
        #     sexActs = self.runRandom.randrange(6,11,1)
        # elif diceroll < 0.019 + 0.082 + 0.063 + 0.072:
        #     sexActs = self.runRandom.randrange(12,23,1)
        # elif diceroll < 0.019 + 0.082 + 0.063 + 0.072 + 0.299:
        #     sexActs = self.runRandom.randrange(24,35,1)
        # elif diceroll < 0.019 + 0.082 + 0.063 + 0.072 + 0.299 + 0.200:
        #     sexActs = self.runRandom.randrange(36,51,1)
        # elif diceroll < 0.019 + 0.082 + 0.063 + 0.072 + 0.299 + 0.200 + 0.124:
        #     sexActs = self.runRandom.randrange(52,155,1)
        # else: #0.141
        #     sexActs = 156
        #
        # return sexActs


    #@profile
    def _become_HIV(self, agent, time):
        """
        :Purpose:
            agent becomes HIV agent. Update all appropriate list and
            dictionaries.

        :Input:
            agent : int

        """
        #print "\t!\tAgent %d just became HIV" % agent
        #agent_cl = self.All_agentSet.get_agent(agent)
        if agent._HIV_bool: pass
        else:
            agent._HIV_bool = True
            agent._HIV_time = 1
            self.NewInfections.add_agent(agent)
            # print "\t\t\t\tAgent %d added to new infection list"%agent.get_ID()
            self.HIV_agentSet.add_agent(agent)
            if agent._PrEP_time > 0:
                if self.runRandom.random() < params.PrEP_resist:
                    agent._PrEPresistance = 1

        if agent._PrEP_bool:
            self._discont_PrEP(agent, time, force=True)

        #self.All_agentSet._subset["HIV"].add_agent(agent)
        #print "\t\t\t\t\tNew HIV agent %d     VL:%.1lf    time:%d"%(agent, Viral_load, time)


    def _sex_possible(self, agent_sex_type, partner_sex_type):
        """
        :Purpose:
        Determine if sex is possible.

        :Input:
        agent_sex_type : str

        partner_sex_type : str

        :Output:
        SexPossible : bool
        """

        # Check input
        if agent_sex_type not in params.agentSexTypes:
            raise ValueError("Invalid agent_sex_type! %s" % str(agent_sex_type))
        if partner_sex_type not in params.agentSexTypes:
            raise ValueError("Invalid partner_sex_type! %s" % str(
                partner_sex_type))

        # Sex possible
        if agent_sex_type == 'HM' and partner_sex_type in ['HF', 'WSW', 'MTF']:
            SexPossible = True
        #elif partner_sex_type == 'HM' and agent_sex_type in ['HF', 'WSW']:
        #    SexPossible = True
        elif agent_sex_type == 'MSM' and partner_sex_type in ['MSM', 'WSW', 'HF', 'MTF']:
            SexPossible = True
        #elif partner_sex_type == 'MSM' and agent_sex_type in ['MSM', 'WSW', 'HF']:
        #    SexPossible = True
        elif agent_sex_type == 'WSW' and partner_sex_type in ['MSM', 'WSW', 'HM']:
            SexPossible = True
        #elif partner_sex_type == 'WSW' and agent_sex_type in ['MSM', 'WSW', 'HM']:
        #    SexPossible = True
        elif agent_sex_type == 'HF' and partner_sex_type in ['HM', 'MSM']:
            SexPossible = True
        elif agent_sex_type == 'MTF' and partner_sex_type in ['HM', 'MSM']:
            SexPossible = True
        else:
            SexPossible = False

        if agent_sex_type == 'HM' and partner_sex_type == 'HM' and SexPossible:
            raise ValueError("Check _sex_possible method!")

        return SexPossible


    def _drug_cessation(self, agent, agent_drug_type):
        """
        :Purpose:
            Account for drug cessation of IDU to NIDU and NIDU to ND.

        :Input:
            agent : int

        """
        if agent_drug_type == 'IDU':
            if self.runRandom.random() < 0.017 / 12 / 2:
                self.tmp_Agents[agent].update({'Drug Type': 'NIDU'})
                if agent in self.tmp_IDU_agents:  # agent might have transitioned already
                    self.tmp_IDU_agents.remove(agent)
                if agent not in self.tmp_NIDU_agents:
                    self.tmp_NIDU_agents.append(agent)
        elif agent_drug_type == 'NIDU':
            if self.runRandom.random() < 0.017 / 12:
                self.tmp_Agents[agent].update({'Drug Type': 'ND'})
                if agent in self.tmp_NIDU_agents:  # agent might have transitioned already
                    self.tmp_NIDU_agents.remove(agent)
                if agent not in self.tmp_ND_agents:
                    self.tmp_ND_agents.append(agent)
                if agent in self.DrugTreatmentAgents_current:
                    self._exit_drug_treatment(agent)
        else:
            raise ValueError('Drug cessation only valid for IDU and NIDU!')

    def _enroll_treatment(self, time):
        """
        :Purpose:
            Account for drug cessation of IDU to NIDU and NIDU to ND.

        :Input:
            agent : int

        """
        print("\n\n!!!!Engaginge treatment process: %d"%time)
        self.treatmentEnrolled = True
        for agent in self.All_agentSet.iter_agents():#self.Agents: #NEW METHOD
            #agent.print_agent()
            #print("\nAgent:",agent)
            # agent_dict = self.Agents[agent]
            if self.runRandom.random() < params.treatmentCov and agent._DU == 'IDU':
                agent._SNE_bool = True




    def _incarcerate(self, agent, time):
        """
        :Purpose:
            Account for drug cessation of IDU to NIDU and NIDU to ND.

        :Input:
            agent : int

        """

        #agent_dict = self.Agents[agent]
        drug_type = agent._DU#self.get_agent_characteristic(agent, 'Drug Type')
        sex_type = agent._SO#self.get_agent_characteristic(agent, 'Sex Type')
        race_type = agent._race#self.get_agent_characteristic(agent, 'Race')
        hiv_bool = agent._HIV_bool#self.get_agent_characteristic(agent, 'HIV')
        tested = agent._tested#self.get_agent_characteristic(agent, 'Tested')
        incar_t = agent._incar_time#agent_dict['incar_t']
        incar_bool = agent._incar_bool
        haart_bool = agent._HAART_bool

        if incar_bool:#agent in self.Incarcerated:
            agent._incar_time -= 1


            #get out if t=0
            if incar_t == 1: #FREE AGENT
                self.incarcerated_agentSet.remove_agent(agent)
                self.NewIncarRelease.add_agent(agent)
                agent._incar_bool = False
                agent._ever_incar_bool = True
                if not agent._highrisk_bool:        #If behavioral treatment on and agent HIV, ignore HR period.
                    if params.inc_treat_HRsex_beh and hiv_bool and (time >=params.inc_treatment_startdate):
                        pass
                    else:                           #Else, become high risk
                        self.highrisk_agentsSet.add_agent(agent)
                        if not agent._everhighrisk_bool:
                            self.NewHRrolls.add_agent(agent)

                        agent._mean_num_partners = agent._mean_num_partners + params.HR_partnerScale
                        agent._highrisk_bool = True
                        agent._everhighrisk_bool = True
                        agent._highrisk_time = params.HR_M_dur


                if (params.inc_treat_RIC or params.inc_treat_HRsex_beh or params.inc_treat_IDU_beh) and (time >=params.inc_treatment_startdate):
                    agent._incar_treatment_time = params.inc_treatment_dur

                if hiv_bool:
                    if haart_bool:
                        if self.runRandom.random() > params.inc_ARTdisc: #12% remain surpressed
                            pass

                        else:
                            agent._HAART_bool = False
                            agent._HAART_adh = 0
                            self.Trt_ART_agentSet.remove_agent(agent)

                        ### END FORCE ####

        elif self.runRandom.random() < params.DemographicParams[race_type][sex_type]['INCAR'] * (1+(hiv_bool*4)) * params.cal_IncarP:
            toss = 2#random.choice( (1, 2) )
            if toss == 1: #JAIL
                timestay = self.runRandom.randint(params.inc_JailMin, params.inc_JailMax) #int(random.triangular(8, 21, 15))
                if hiv_bool and not tested:
                    if self.runRandom.random() < params.inc_JailTestProb :
                        agent._tested = True #self.tmp_Agents[agent].update({'Tested': 1})
                        #self.HIVidentified_agents.append(agent)

            else: #PRISON
                timestay = self.runRandom.randint(params.inc_PrisMin, params.inc_PrisMax)
                if hiv_bool:
                    if not tested:
                        if self.runRandom.random() < params.inc_PrisTestProb:
                            agent._tested = True #self.tmp_Agents[agent].update({'Tested': 1})
                            #self.HIVidentified_agents.append(agent)
                    else: #Then tested and HIV, check to enroll in ART
                        if self.runRandom.random() < params.inc_ARTenroll:
                            tmp_rnd = self.runRandom.random()
                            HAART_ADH = params.inc_ARTadh
                            if tmp_rnd < HAART_ADH:
                                adherence = 5
                            else:
                                adherence = self.runRandom.randint(1,4)

                            #Add agent to HAART class set, update agent params
                            agent._HAART_bool = True
                            agent._HAART_adh = adherence
                            agent._HAART_time = time
                            self.Trt_ART_agentSet.add_agent(agent)

            agent._incar_bool = True
            agent._incar_time = timestay
            self.incarcerated_agentSet.add_agent(agent)
            self.totalIncarcerated += 1

            #PUT PARTNERS IN HIGH RISK
            for tmpA in agent._partners:
                if tmpA._highrisk_bool == True:
                    pass
                    #print "ALREADY HR"
                else:
                    if self.runRandom.random() < params.HR_proportion:
                        #print "Making agent %d (%s) HR"%(tmpA._ID, tmpA._SO)
                        if not tmpA._highrisk_bool:
                            self.highrisk_agentsSet.add_agent(tmpA)
                            if not tmpA._everhighrisk_bool:
                                self.NewHRrolls.add_agent(tmpA)
                            tmpA._mean_num_partners += params.HR_partnerScale #32.5 #2 + 3.25 from incar HR
                            tmpA._highrisk_bool = True
                            tmpA._everhighrisk_bool = True
                            tmpA._highrisk_time = params.HR_F_dur


    def _HIVtest(self, agent, time):
        """
        :Purpose:
            Test the agent for HIV. If detected, add to identified list.

        :Input:
            agent : agent_Class
            time : int

        :Output:
            none
        """
        # Drug Treatment
        # SEP

        drug_type = agent._DU#self.get_agent_characteristic(agent, 'Drug Type')
        sex_type = agent._SO#self.get_agent_characteristic(agent, 'Sex Type')
        race_type = agent._race#self.get_agent_characteristic(agent, 'Race')
        hiv_Status = agent._HIV_bool#self.get_agent_characteristic(agent, 'HIV')
        tested = agent._tested#self.get_agent_characteristic(agent, 'Tested')
        if not tested:
            test_prob = params.DemographicParams[race_type][sex_type]['HIVTEST']

            #Rescale based on calibration param
            test_prob = test_prob *  params.cal_TestFreq

            #If roll less than test probablity
            if self.runRandom.random() < test_prob: ###WAS / 10
                # Become tested, add to tested agent set
                agent._tested = True
                self.NewDiagnosis.add_agent(agent)
                self.Trt_Tstd_agentSet.add_agent(agent)
                # If treatment co-enrollment enabled and coverage greater than 0
                if self.treatmentEnrolled and params.treatmentCov > 0:
                    #For each partner, attempt to test for HIV
                    for ptnr in agent._partners:
                        if ptnr._HIV_bool and not ptnr._tested:
                            if self.runRandom.random < 0.87:
                                ptnr._tested = True
                                self.NewDiagnosis.add_agent(ptnr)





    def _VCT(self, agent, time):
        """
        :Purpose:
            Account for voluntary Counseling and Testing(VCT)

        :Input:
            agent : int
            partner : int
            time : int

        :Output:
            none
        """
        # Drug Treatment
        # SEP
        drug_type = self.get_agent_characteristic(agent, 'Drug Type')
        SEPstat = False
        if agent in self.SEPAgents:
            if time == self.SEPAgents[agent]:
                SEPstat = True
        if drug_type == 'IDU':
            if SEPstat:
                if self.runRandom.random() < self.VCT_NSP:  # !!!!!!!!!!!!!!!!!!!!
                    self.VCTAgents.update({agent: time})
            else:
                if self.runRandom.random() < self.VCT_NoNSP_IDU:  # !!!!!!!!!!!!!!!!!!!
                    self.VCTAgents.update({agent: time})
        if drug_type == 'NIDU':
            if self.runRandom.random() < self.VCT_NoNSP_NIDU:  # !!!!!!!!!!!!!!!!!!
                self.VCTAgents.update({agent: time})
        elif agent in self.MSM_agents and self.runRandom.random() < self.VCT_NoNSP_MSM:  # !
            self.VCTAgents.update({agent: time})
        else:
            if self.runRandom.random() < self.VCT_NoNSP_EE:  # !!!!!!!!!!!!!!!!!!!!!!!!!
                self.VCTAgents.update({agent: time})


    def _SEP(self, agent, time):
        """
        :Purpose:
            Account for SEP (Syringe Exchange Program) for IDU agents. \n

        :Input:
            time : int

        :Output:
            SEPstat : bool
        """
        if self.get_agent_characteristic(agent, 'Drug Type') != 'IDU':
            raise ValueError("_SEP only valid for IDU agents! agent: %d %s" %
                             (agent, str(self.get_agent_characteristic(agent,
                                                                       'Drug Type'))))
        # Drug treatment increases likelihood of SEP use
        if agent in self.DrugTreatmentAgents_current and self.IDU_agents:

            if self.runRandom.random() < self.NSP_SAT:
                self.SEPAgents.update({agent: time})
                SEPstat = True
            else:
                SEPstat = False
        elif self.runRandom.random() < self.NSP_NoSAT:
            self.SEPAgents.update({agent: time})
            SEPstat = True
        else:
            SEPstat = False

        return SEPstat


    def _exit_drug_treatment(self, agent):
        """
        Agent exits drug treament.
        """
        self.DrugTreatmentAgents_past.update({agent:
                                                  self.DrugTreatmentAgents_current[agent]})
        del self.DrugTreatmentAgents_current[agent]


    def _enter_and_exit_drug_treatment(self, agent, time):
        """
        :Purpose:
            Account for drug treatment for IDU agents. \n
            Entering drug treatment is similar to SEP, drug treatment has
            a functional relationship given as follows:
            P(treatment (t+1) | IDU or NIDU) = P(treatment (agent,t) | IDU) if x < N
            else: 0 (x >= 0)
            where N is the total number of treatment slots available. \n
            An agent who was already in drug treatment and relapsed, has a
            pobability twice as strong to reenter drug treatment
            at a later point.

        :Input:
            agent : int

            time : int

        :Output:
            bool
        """
        agent_drug_type = self.get_agent_characteristic(agent, 'Drug Type')

        N_TrSpots_Max = 100000  # max number of treatment spots
        if (agent in self.DrugTreatmentAgents_current and self.runRandom.random() < self.SAT_disc):
            self._exit_drug_treatment(agent)
        elif self.N_TreatmentSpots < N_TrSpots_Max:
            if agent in self.SEPAgents:
                prob = self.SAT_NSP
            elif agent in self.DrugTreatmentAgents_past:
                prob = 0.18
            else:
                if agent_drug_type == 'IDU':
                    prob = self.SAT_NoNSP_IDU
                elif agent_drug_type == 'NIDU':
                    prob = self.SAT_NIDU
                else:
                    mssg = 'Drug treatment only valid for NIDU and IDUs! %s'
                    raise ValueError(mssg % agent_drug_type)
            if self.runRandom.random() < prob and agent not in self.tmp_ND_agents:
                self.DrugTreatmentAgents_current.update({agent: time})
                self.N_TreatmentSpots += 1
        elif self.N_TreatmentSpots == N_TrSpots_Max:
            pass
        else:
            mssg = 'Check self.N_TreatmentSpots! Max value = 500! %d'
            raise ValueError(mssg % self.N_TreatmentSpots)


    def _initiate_HAART(self, agent, time):
        """
        :Purpose:
            Account for HIV treatment through highly active antiretroviral therapy (HAART).
            HAART was implemented in 1996, hence, there is treatment only after 1996.
            HIV treatment assumes that the agent knows his HIV+ status.

        :Input:
            time : int

        :Output:
            none
        """

        HAART_coverage = 0.30

        # Check valid input
        if not agent._HIV_bool:
            print "HIV_agents: ", sorted(self.HIV_agents)
            print "tmp_HIV_agents: ", sorted(self.tmp_HIV_agents)
            print "Agent[agent]", self.Agents[agent]
            try:
                print "tmp_Agent[agent]", self.tmp_Agents[agent]
            except KeyError:
                pass
            raise ValueError("HAART only valid for HIV agents!agent:%s" %
                             str(agent))

        agent_drug_type = agent._DU#self.get_agent_characteristic(agent, 'Drug Type')
        agent_haart = agent._HAART_bool#self.get_agent_characteristic(agent, 'HAARTa')
        agent_HIV = agent._HIV_bool#self.get_agent_characteristic(agent, 'HIV')
        agent_Test_bool = agent._tested#self.Agents[agent]['Tested']
        agent_race = agent._race
        agent_so = agent._SO

        #Set HAART agents adherence at t=0 for all instanced HAART
        if time == 0 and agent_haart:
            #agent_haart = agent._HAART_bool#self.get_agent_characteristic(agent, 'HAARTa')
            agent_haart_adh = agent._HAART_adh
            if agent_haart_adh == 0:
                tmp_rnd = self.runRandom.random()
                HAART_ADH = params.DemographicParams[agent_race][agent_so]['HAARTadh']
                if tmp_rnd < HAART_ADH:
                    adherence = 5
                else:
                    adherence = self.runRandom.randint(1,4)

                #add to agent haart set
                agent._HAART_adh = adherence
                agent._HAART_time = time
                #self.Trt_ART_agentSet.add_agent(agent)


        # Determine probability of HIV treatment
        if time >= 0 and agent_Test_bool:
            #OLD WAY OF PUTTING AGENTS ON HAART
            prob = 0
            if agent_Test_bool:
                if agent_drug_type == 'IDU':
                    if agent:# in self.DrugTreatmentAgents_current:
                        prob = 0.00625

                elif agent_drug_type == 'NIDU':
                    if agent in self.DrugTreatmentAgents_current:
                        prob = 0.0117

                elif agent_drug_type == 'NDU':
                    prob = 0.0117
            else:
                prob = 0.0

            # Go on HAART
            if not agent_haart and agent._HAART_time == 0:
                if self.runRandom.random() < prob * params.cal_ART_cov:

                    #self.tmp_Agents[agent].update({'HAARTa': 1})
                    #self.tmp_HAART_agents.append(agent)

                    tmp_rnd = self.runRandom.random()
                    HAART_ADH = params.DemographicParams[agent_race][agent_so]['HAARTadh']
                    if tmp_rnd < HAART_ADH:
                        adherence = 5
                    else:
                        adherence = self.runRandom.randint(1,4)

                    #Add agent to HAART class set, update agent params
                    agent._HAART_bool = True
                    agent._HAART_adh = adherence
                    agent._HAART_time = time
                    self.Trt_ART_agentSet.add_agent(agent)

                    #self.AdherenceAgents.update({agent: adherence})

            elif agent_haart and self.runRandom.random() < params.DemographicParams[agent_race][agent_so]['HAARTdisc']:
                if agent._incar_treatment_time > 0 and params.inc_treat_RIC:
                    pass
                else:
                    agent._HAART_bool = False
                    agent._HAART_adh = 0
                    agent._HAART_time = 0
                    self.Trt_ART_agentSet.remove_agent(agent)



    def _PrEP_elligible(self, agent, time):
        elligble = False
        if params.PrEP_target_model == 'Allcomers':
            elligble = True
        elif params.PrEP_target_model == 'HighPN5':
            if agent._mean_num_partners >= 5:
                elligble = True
        elif params.PrEP_target_model == 'HighPN10':
            if agent._mean_num_partners >= 10:
                elligble = True
        elif params.PrEP_target_model == 'SRIns':
            if agent._sexualRole == 'Insertive':
                elligble = True
        elif params.PrEP_target_model == 'MSM':
            if agent._SO == ('MSM' or 'MTF'):
                elligble = True
        elif params.PrEP_target_model == 'RandomTrial':
            # If using random trial
            if time == 0:
                #if in init timestep 0, use agent set elligiblity
                elligible = agent._PrEP_elligible
            if time > 0:
                #else, false to not continue enrollment past random trial start
                elligble = False

        return elligble


    def _calc_PrEP_load(self, agent):
        """
        :Purpose:
            Determine load of PrEP concentration in agent.

        :Input:
            agent : agent()

        :Output:
            none
        """

        #N(t) = N0 (0.5)^(t/t_half)
        #print agent._PrEP_load, agent._PrEP_lastDose

        agent._PrEP_lastDose += 1
        if agent._PrEP_lastDose > 12:
            agent._PrEP_load = 0.0
        else:
            agent._PrEP_load = params.PrEP_peakLoad * ((0.5)**(agent._PrEP_lastDose/(params.PrEP_halflife/30)))
        #print agent._ID, agent._PrEP_load, agent._PrEP_lastDose


    def _discont_PrEP(self, agent, time, force=False):

        #If force flag set, auto kick off prep.
        if force==True:
            self.Trt_PrEP_agentSet.remove_agent(agent)
            agent._PrEP_bool = False
        #else if agent is no longer enrolled on PrEP, increase time since last dose
        elif agent._PrEP_time > 0:
            if agent._PrEP_time == 1:
                agent._PrEP_bool = False
                agent._PrEP_time -= 1
            else:
                agent._PrEP_time -= 1

        #else if agent is on PrEP, see if they should discontinue
        elif agent._PrEP_bool and agent._PrEP_time == 0:
            if self.runRandom.random() < params.DemographicParams[agent._race][agent._SO]['PrEPdisc']:
                agent._PrEP_time = params.PrEP_falloutT
                self.Trt_PrEP_agentSet.remove_agent(agent)

                if params.PrEP_type == 'Oral':
                    agent._PrEP_bool = False
                #print 'Agent%d removed from PrEP\tP_T:%d'%(agent._ID, agent._PrEP_time)
            else: #if not discontinue, see if its time for a new shot.
                if agent._PrEP_lastDose > 2:
                    #print 'Agent%d renewed PrEP'%(agent._ID)
                    agent._PrEP_lastDose = -1

        if params.PrEP_type == 'Inj':
            self._calc_PrEP_load(agent)

    def _initiate_PrEP(self, agent, time, force=False):
        """
        :Purpose:
            Place agents onto PrEP treatment.
            PrEP treatment assumes that the agent knows his HIV+ status is negative.

        :Input:
            time : int

        :Output:
            none
        """
        def _enrollPrEP(self, agent):
            agent._PrEP_bool = True
            agent._PrEP_time = 0
            self.Trt_PrEP_agentSet.add_agent(agent)
            tmp_rnd = self.runRandom.random()

            if tmp_rnd < params.PrEP_Adherence:
                agent._PrEP_adh = 1
            else:
                agent._PrEP_adh = 0

            #set PrEP load and dosestep for PCK
            if params.PrEP_type == 'Inj':
                agent._PrEP_load = params.PrEP_peakLoad
                agent._PrEP_lastDose = 0


        if agent == None:
            print "OHHH boi no prep agent"
            return None
        # Check valid input
        #agent.print_agent()
        if agent._PrEP_bool:
            print agent._PrEP_bool
            raise ValueError("PrEP only valid for agents not on PrEP!agent:%d" % agent.get_ID())

        if agent._HIV_bool:
            raise ValueError("PrEP only valid for HIV- agents!agent:%d" % agent.get_ID())

        # Determine probability of HIV treatment
        agent_drug_type = agent._DU
        agent_race = agent._race
        agent_so = agent._SO

        if force:
            _enrollPrEP(self, agent)
        else:
            numPrEP_agents = self.Trt_PrEP_agentSet.num_members()

            if params.PrEP_target_model == 'Clinical':
                target_PrEP_population = self.All_agentSet.num_members() - self.HIV_agentSet.num_members()
                target_PrEP = target_PrEP_population * params.PrEP_Target
            else:
                target_PrEP = int((self.All_agentSet.num_members()-self.All_agentSet._subset["HIV"].num_members()) * params.PrEP_Target)

            if params.PrEP_clinic_cat = 'Racial' and agent_race == 'BLACK':
                if self.runRandom.random() < params.PrEP_Target:
                    _enrollPrEP(self, agent)
            elif numPrEP_agents < target_PrEP and time >= params.PrEP_startT:
                #print 'Agent%d added from PrEP'%(agent._ID)
                _enrollPrEP(self, agent)
            

    def _get_clinic_agent(self, clinicBin, elligiblePool):
        i=1
        pMatch = params.clinicAgents[clinicBin][i]['Prob']
        RN = self.runRandom.random()
        while(True):
            if RN <= pMatch:
                break
            else:
                i+=1
                pMatch += params.clinicAgents[clinicBin][i]['Prob']
            if i==5:break

        minNum = params.clinicAgents[clinicBin][i]['min']
        maxNum = params.clinicAgents[clinicBin][i]['max']

        iterations = 1
        while iterations < 3:
            randomK_sample = self.runRandom.sample(elligiblePool,params.cal_ptnrSampleDepth)
            #randomK_sample = self.runRandom.sample(self.All_agentSet._subset["MSM"]._members,params.cal_ptnrSampleDepth)
            elligibleK_Pool = [ag for ag in randomK_sample if ((ag._mean_num_partners >= minNum) and (ag._mean_num_partners <= maxNum))]
            #for a in elligibleK_Pool: print minNum, maxNum, a._mean_num_partners
            if elligibleK_Pool:
                selected = self.runRandom.choice(elligibleK_Pool)
                elligiblePool.remove(selected)
                return selected
            else:
                print "Looking for agent with min:%d and max %d failed %d times"%(minNum, maxNum,iterations)
                iterations += 1
                #print "Looking in another subsetK"

        print "No suitable PrEP agent"
        #raise ValueError("No suitable PrEP agent")
        return None



    def _progress_to_AIDS(self, agent, agent_drug_type):
        """
        :Purpose:
            Model the progression of HIV agents to AIDS agents
        """
        # only valid for HIV agents
        if not agent._HIV_bool:
            raise ValueError("HAART only valid for HIV agents!agent:%s" % str(agent._ID))

        #if agent not in self.AIDS_agents:
        if not agent._HAART_bool:
            adherenceStat = agent._HAART_adh#self.AdherenceAgents[agent]
            if adherenceStat > 0:
                if adherenceStat == 1:
                    prob = 0.0051

                if adherenceStat == 2:
                    prob = 0.0039

                if adherenceStat == 3:
                      prob = 0.0032

                if adherenceStat == 4:
                    prob = 0.0025

                if adherenceStat == 5:
                    prob = 0.0008

            else:
                prob = 0.0051

            if self.runRandom.random() < prob * params.cal_ProgAIDS:
                #print "--------AIDS? -> yes"
                agent._AIDS_bool = True
                self.HIV_AIDS_agentSet.add_agent(agent)
                #self.tmp_AIDS_agents.append(agent)
                #self.tmp_Agents[agent].update({'AIDS': 1})


    def _reset_death_count(self):
        self.num_Deaths = {}
        for HIV_status in ['Total','HIV-', 'HIV+']:
            self.num_Deaths.update({HIV_status: {}})
            for tmp_type in [HIV_status, 'MSM', 'HM', 'HF', 'WSW', 'MTF']:
                self.num_Deaths[HIV_status].update({tmp_type: 0})


    def _remove_agent(self, agent):
        """
        :Purpose:
            Remove agent from the population.
            Delete agent is a key to an associated dictionary which stores the internal.

        :Input:
            agent : int

        """

        # Drug type lists (agent might have been updated)
        drug_type = self.tmp_Agents[agent]['Drug Type']

        if drug_type == 'IDU':
            try:
                self.tmp_IDU_agents.remove(agent)
            except ValueError:
                pass
        elif drug_type == 'NIDU':
            try:
                self.tmp_NIDU_agents.remove(agent)
            except ValueError:
                pass
        elif drug_type == 'ND':
            try:
                self.tmp_ND_agents.remove(agent)
            except ValueError:
                pass
                #print "Agents[agent]", self.Agents[agent]
                #print "tmp_Agents[agent]", self.tmp_Agents[agent]
                #print "tmp_ND_agents", sorted(self.tmp_ND_agents)
        else:
            raise ValueError("Invalid drug type! %s" % str(drug_type))

        # Sex type lists
        try:
            # Agent might have been updated
            sex_type = self.tmp_Agents[agent]['Sex Type']
        except KeyError:
            sex_type = self.get_agent_characteristic(agent, 'Sex Type')
        if sex_type == 'MSM':
            self.tmp_MSM_agents.remove(agent)
        elif sex_type == 'HF':
            self.tmp_HF_agents.remove(agent)
        elif sex_type == 'WSW':
            self.tmp_WSW_agents.remove(agent)
        elif sex_type == 'HM':
            self.tmp_HM_agents.remove(agent)
        else:
            raise ValueError("Invalid sex type! %s" % str(sex_type))

        # HIV and AIDS lists
        try:
            self.tmp_HIV_agents.remove(agent)
        except ValueError:
            pass

        try:
            self.tmp_AIDS_agents.remove(agent)
        except ValueError:
            pass

        try:
            self.Incarcerated.remove(agent)
        except ValueError:
            pass

        try:
            self.tmp_HAART_agents.remove(agent)
        except ValueError:
            #print "WTFFFFF"
            pass

        try:
            self.HIVidentified_agents.remove(agent)
        except ValueError:
            #print "WTFFFFF"
            pass

        # Other lists / dictionaries
        if agent in self.SEPAgents:
            # dict of users who used SEP (agent:time)
            del self.SEPAgents[agent]
        if agent in self.SEPAgents_past:
            del self.SEPAgents_past[agent]
        if agent in self.DrugTreatmentAgents_current:
            # dictionary of users who are currently undergoing
            del self.DrugTreatmentAgents_current[agent]
        if agent in self.DrugTreatmentAgents_past:
            # dictionary of users who underwent drug treatment
            del self.DrugTreatmentAgents_past[agent]
        if agent in self.VCTAgents:
            # list of agents who get tested for HIV ((agent:time)
            del self.VCTAgents[agent]
        for time in self.HIV_key_transitiontime:
            tmp_agents = self.HIV_key_transitiontime[time]
            if agent in tmp_agents:
                tmp_agents.remove(agent)
            self.HIV_key_transitiontime.update({time: tmp_agents})


    def _die_and_replace(self):
        """
        :Purpose:
            Let agents die and replace the dead agent with a new agent randomly.
        """
        totalDeaths = 0
        # self.num_Deaths["Total"] = 0
        # self.num_Deaths["HIV+"] = 0
        # self.num_Deaths["HIV-"] = 0
        #dynnetworkReport = open('Results/dynnetworkReport.txt', 'a')
        for agent in self.All_agentSet._members:#iter_agents(): #self.Agents:

            if agent._incar_bool:#self.IncarceratedClass.is_member(agent):
                #print "Agent %d is incarcerated. Cannot die" % agent.get_ID()
                pass
            else:
                # Probability for dying
                drug_type = agent._DU #self.get_agent_characteristic(agent, 'Drug Type')
                sex_type = agent._SO #self.get_agent_characteristic(agent, 'Sex Type')
                HIV_status = agent._HIV_bool #self.get_agent_characteristic(agent, 'HIV')
                AIDSStatus = agent._AIDS_bool #self.get_agent_characteristic(agent, 'HIV')
                agent_Race = agent._race #self.get_agent_characteristic(agent, 'Race')
                #adherence =  self.AdherenceAgents[agent]


                if HIV_status:
                    if AIDSStatus: #AIDS DEATH RATE
                        if agent_Race == 'WHITE':
                            p = 34.4
                        elif agent_Race == 'BLACK':
                            p = 41.6
                        else:
                            raise ValueError("Invalid RACE type! %s" % str(agent_Race))
                        #p = self.ProbDeath[drug_type]['AIDS']

                    elif agent._HAART_adh > 1: #HAART DEATH RATE
                        if agent_Race == 'WHITE':
                            p = 8.6
                        elif agent_Race == 'BLACK':
                            p = 10.4
                        else:
                            raise ValueError("Invalid RACE type! %s" % str(agent_Race))
                        #p = self.ProbDeath[drug_type]['HIV+/HAART']

                    else: #HIV+ DEATH RATE
                        if agent_Race == 'WHITE':
                            p = 17.2
                        elif agent_Race == 'BLACK':
                            p = 20.8
                        else:
                            raise ValueError("Invalid RACE type! %s" % str(agent_Race))
                        #p = self.ProbDeath[drug_type]['HIV+']
                    p = p * params.cal_Mortality

                elif not HIV_status: # NON HIV DEATH RATE
                    if agent_Race == 'WHITE':
                        p = 8.6
                    elif agent_Race == 'BLACK':
                        p = 10.4
                    else:
                        raise ValueError("Invalid RACE type! %s" % str(agent_Race))
                    #p = self.ProbDeath[drug_type]['HIV-']

                else:
                    raise ValueError("Invalid HIV type! %s" % str(HIV_status))

                #print p
                p = p / 12000.0#12000.0 #putting it into per 1 person-month
                if self.runRandom.random() < p:
                    # print "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tAgent %d died rolling under %.10lf" % (agent.get_ID(), p)
                    
                    totalDeaths += 1
                    if HIV_status: ident = "HIV+"
                    else: ident = "HIV-"
                    self.num_Deaths["Total"][sex_type] += 1
                    self.num_Deaths[ident][sex_type] += 1
                    ID_number = agent.get_ID()
                    race = agent._race

                    #End all existing relationships
                    for rel in agent._relationships:
                        # print "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDeleting relationship between %d and %d" % (rel._ID1.get_ID(), rel._ID2.get_ID())
                        rel.progress(forceKill=True)

                        self.Relationships.remove_agent(rel)

                    
                    #Remove agent node and edges from network graph
                    self.get_Graph().remove_node(agent)



                    #Remove agent from agent class and sub-sets
                    self.All_agentSet.remove_agent(agent)

                    #Delete agent object
                    del agent

                    # Create new agent
                    agent_cl = self._return_new_Agent_class(ID_number,race)
                    self.create_agent(agent_cl, race)
                    # print("Adding node {}".format(agent_cl))

                    self.G.add_node(agent_cl)

                elif 1==0:
                    # Replace with new agent (random characteristics)
                    rv = self.runRandom.random()
                    if rv < params.DemographicParams['WHITE']['ALL']['Proportion']:
                        deliminator = 'WHITE'
                        #print "\t\tReplaced with ND"
                    else:
                        deliminator = 'BLACK'
                        #print "\t\tReplaced with NIDU"

                        """#################### OLD WAY
                    # Replace with new agent (random characteristics)
                    rv = self.runRandom.random()
                    if rv < 0.9229:
                        drug_type = 'ND'
                        print "\t\tReplaced with ND"
                    elif rv < 0.9229 + 0.0647:
                        drug_type = 'NIDU'
                        print "\t\tReplaced with NIDU"
                    else:
                        drug_type = 'IDU'
                        print "\t\tReplaced with "
                        """


                     #################### NOW SET TO REPLAC WITH WHAT DIED"

                    #drug_type = 'IDU'

                    # New agent dict
                    agent_dict = self._return_new_agent_dict(deliminator)
                    if deliminator != agent_dict['Race']:
                        raise ValueError("Inconsistent drug type!%s" % str(agent_dict['Drug Type']))

                    # Update tmp_Agents dictionary with new agent
                    self.tmp_Agents.update({agent: agent_dict})
                    #print "New agent updated"

                    # Drug Type
                    drug_type = agent_dict['Drug Type']
                    if drug_type == 'IDU':
                        self.tmp_IDU_agents.append(agent)
                    elif drug_type == 'NIDU':
                        self.tmp_NIDU_agents.append(agent)
                    elif drug_type == 'ND':
                        self.tmp_ND_agents.append(agent)
                    else:
                        raise ValueError("Invalid drug type! %s" % str(drug_type))

                    # Sex Type
                    SexType = agent_dict['Sex Type']
                    if SexType == 'HM':
                        self.tmp_HM_agents.append(agent)
                    elif SexType == 'HF':
                        self.tmp_HF_agents.append(agent)
                    elif SexType == 'MSM':
                        self.tmp_MSM_agents.append(agent)
                    elif SexType == 'WSW':
                        self.tmp_WSW_agents.append(agent)
                    else:
                        raise ValueError("Invalid SexType! %s" % str(SexType))

                    # HIV
                    HIVStatus = agent_dict['HIV']
                    if HIVStatus == 1:
                        #print "NEW AGENT %d %s WAS HIV"%(agent, drug_type)
                        self.tmp_HIV_agents.append(agent)
                    elif HIVStatus != 0:
                        raise ValueError("Invalid HIVType! %s" % str(HIVStatus))
                    #else:
                        #print "NEW AGENT %d %s WAS NOT HIV"%(agent, drug_type)

                    # AIDS
                    AIDSStatus = agent_dict['AIDS']
                    if AIDSStatus == 1:
                        #print "NEW AGENT WAS AIDS"
                        self.tmp_AIDS_agents.append(agent)
                    elif AIDSStatus != 0:
                        raise ValueError("Invalid AIDS Status! %s" % str(AIDSStatus))

                    # HAART
                    HAARTStatus = agent_dict['HAARTa']
                    if HAARTStatus == 1:
                        #print "NEW AGENT WAS HAART"
                        self.tmp_HAART_agents.append(agent)
                    elif HAARTStatus != 0:
                        raise ValueError("Invalid HAART Status! %s" % str(HAARTStatus))

                    #Incarcerated
                    IncarceratedTime = agent_dict['incar_t']
                    if IncarceratedTime >= 1:
                        self.Incarcerated.append(agent)
                    elif IncarceratedTime < 0:
                        raise ValueError("Invalid AIDS Status! %s"%str(IncarceratedTime))

                    # Check
                    if HIVStatus == 1:
                        if agent not in self.tmp_HIV_agents:
                            raise ValueError("Agent must be in HIV_agents")
                    if AIDSStatus == 1:
                        if agent not in self.tmp_AIDS_agents:
                            raise ValueError("Agent must be in AIDS_agents")
                    if HAARTStatus == 1:
                        if agent not in self.tmp_HAART_agents:
                            raise ValueError("Agent must be in HAART_agents")

                    """# write new agent to dynnetworkReport
                    #print "Writing to dynNetReport"
                    reportLine = '\t'.join([repr(self.TimeStep), 'DEATH', repr(agent)])
                    #open('Results/dynnetworkReport.txt', 'a').write('\n' + reportLine)
                    dynnetworkReport.write('\n' + reportLine)

                    reportLine = '\t'.join([repr(self.TimeStep), 'NEWAGENT', repr(agent), SexType, drug_type, repr(HIVStatus)])
                    #open('Results/dynnetworkReport.txt', 'a').write('\n' + reportLine)
                    dynnetworkReport.write('\n' + reportLine)"""

            #if len(self.Agents) != self.PopulationSize:
            #    raise ValueError("Wrong Population size!%s != %s" % (str(len(self.Agents)), str(self.PopulationSize)))

            #iduPrec = self.num_Deaths['IDU']['IDU']/float(totalDeaths) if totalDeaths else 0
            #niduPrec = self.num_Deaths['NIDU']['NIDU']/float(totalDeaths) if totalDeaths else 0
            #ndPrec = self.num_Deaths['ND']['ND']/float(totalDeaths) if totalDeaths else 0
            """
            print "\n\t=== DEATH REPORT ==="
            print "\tType\tCount\t%"
            print "\tIDU \t%d\t\t%.5lf" % (self.num_Deaths['IDU']['IDU'], iduPrec)
            print "\tNIDU\t%d\t\t%.5lf" % (self.num_Deaths['NIDU']['NIDU'], niduPrec)
            print "\tND  \t%d\t\t%.5lf" % (self.num_Deaths['ND']['ND'], ndPrec)
            print "\tTotal\t%d\n" % totalDeaths
            """
            #dynnetworkReport.close()

    #@profile
    def _update_population(self):
        """
        :Purpose:
            Update the population. Changes resulting from parsing
            through the agents and applying the update rules are stored
            in :py:attr:`tmp_agent_dict`. This method updates the whole
            population, i.e., it copies changes from the :py:attr:`tmp_agent_dict`
            dictionary and copies it into the :py:attr:`Agents` dictionary.
        :Input:
            none

        :Output:
            none
        """

        # self.Agents = deepcopy(self.tmp_Agents)
        #
        # self.IDU_agents = copy(list(set(self.tmp_IDU_agents)))
        # self.NIDU_agents = copy(list(set(self.tmp_NIDU_agents)))
        # self.ND_agents = copy(list(set(self.tmp_ND_agents)))
        #
        # self.HM_agents = copy(list(set(self.tmp_HM_agents)))
        # self.HF_agents = copy(list(set(self.tmp_HF_agents)))
        # self.MSM_agents = copy(list(set(self.tmp_MSM_agents)))
        # self.WSW_agents = copy(list(set(self.tmp_WSW_agents)))
        #
        # self.AIDS_agents = copy(list(set(self.tmp_AIDS_agents)))
        # self.HIV_agents = copy(list(set(self.tmp_HIV_agents)))
        # self.HAART_agents = copy(list(set(self.tmp_HAART_agents)))

        #print "HIV population:", self.HIV_agents
        self.SEPAgents = {}  # SEP has no memory


    def _check_population(self):
        """
        :Purpose:
            Check consistency of population.
            Only called in unittest.

        """

        # Check consistency of last partners
        if (not (np.all(self.AdjMat.sum(0) == self.AdjMat.conj().sum(0)) and
                     np.all(self.AdjMat.sum(1) == self.AdjMat.conj().sum(1)))):
            raise ValueError("Adjacency matrix not symmetric!")

        # Check consistency of real population
        count_HF = 0
        count_HM = 0
        count_MSM = 0
        count_WSW = 0
        count_ND = 0
        count_NIDU = 0
        count_IDU = 0
        count_HIV = 0
        count_AIDS = 0
        for (agent, d) in self.Agents.iteritems():
            agent_dict = d
            # Sex type
            sex_type = agent_dict['Sex Type']
            if sex_type == 'HF':
                if agent not in self.HF_agents:
                    print self.Agents[agent]
                    raise ValueError("Check agents HF Sex type %d" % agent)
                else:
                    count_HF += 1
            elif sex_type == 'HM':
                if agent not in self.HM_agents:
                    print self.Agents[agent]
                    raise ValueError("Check agents HM Sex type %d" % agent)
                else:
                    count_HM += 1
            elif sex_type == 'MSM':
                if agent not in self.MSM_agents:
                    raise ValueError("Check agents MSM Sex type %d" % agent)
                else:
                    count_MSM += 1
            elif sex_type == 'WSW':
                if agent not in self.WSW_agents:
                    print self.Agents[agent]
                    raise ValueError("Check agents WSW Sex type %d" % agent)
                else:
                    count_WSW += 1
            else:
                raise ValueError("Invalid sex type %s" % str(sex_type))

            # Drug type
            drug_type = agent_dict['Drug Type']
            if drug_type == 'ND':
                if agent not in self.ND_agents:
                    print self.Agents[agent]
                    raise ValueError("Check agents ND Drug type %d" % agent)
                else:
                    count_ND += 1
            elif drug_type == 'NIDU':
                if agent not in self.NIDU_agents:
                    print self.Agents[agent]
                    raise ValueError("Check agents NIDU Drug type %d" % agent)
                else:
                    count_NIDU += 1
            elif drug_type == 'IDU':
                if agent not in self.IDU_agents:
                    print self.Agents[agent]
                    raise ValueError("Check agents IDU Drug type %d" % agent)
                else:
                    count_IDU += 1
            else:
                raise ValueError("Invalid drug type %s" % str(drug_type))

            # HIV
            HIVstatus = agent_dict['HIV']
            if HIVstatus != 0:
                if agent not in self.HIV_agents:
                    print self.Agents[agent]
                    raise ValueError("Check agent HIV %d" % agent)
                else:
                    count_HIV += 1
            # AIDS
            AIDSstatus = agent_dict['AIDS']
            if AIDSstatus != 0:
                if agent not in self.AIDS_agents:
                    print self.Agents[agent]
                    raise ValueError("Check agent AIDS %d" % agent)
                else:
                    count_AIDS += 1

        if len(self.HF_agents) != count_HF:
            raise ValueError("self.HF agents contains too many agents!")
        if len(self.HM_agents) != count_HM:
            print "len(self.HM_agents)=%d" % len(self.HM_agents)
            print "count_HM=%d" % count_HM
            raise ValueError("self.HM agents contains too many agents!")
        if len(self.MSM_agents) != count_MSM:
            raise ValueError("self.MSM agents contains too many agents!")
        if len(self.WSW_agents) != count_WSW:
            raise ValueError("self.WSW agents contains too many agents!")

        if len(self.NIDU_agents) != count_NIDU:
            raise ValueError("self.NIDU_agents contains too many agents!")
        if len(self.ND_agents) != count_ND:
            raise ValueError("self.ND agents contains too many agents!")
        if len(self.IDU_agents) != count_IDU:
            mssg = "self.IDU agents contains too many agents!\
                    \nlen(self.IDU_agents)=%d\ncount_IDU=%d\n"
            raise ValueError(mssg % (len(self.IDU_agents), count_IDU))

        if len(self.HIV_agents) != count_HIV:
            raise ValueError("self.HIV_agents contains too many agents!\
                \nlen(self.HIV_agents) = %d\ncount_HIV = %d\n" % (
                len(self.HIV_agents), count_HIV))
        if len(self.AIDS_agents) != count_AIDS:
            raise ValueError("self.AIDS agents contains too many agents!")

        # Check consistency of tmp population
        count_HF = 0
        count_HM = 0
        count_MSM = 0
        count_WSW = 0
        count_ND = 0
        count_NIDU = 0
        count_IDU = 0
        count_HIV = 0
        count_AIDS = 0
        for (agent, d) in self.tmp_Agents.iteritems():
            agent_dict = d
            # Sex type
            sex_type = agent_dict['Sex Type']
            if sex_type == 'HF':
                if agent not in self.tmp_HF_agents:
                    print self.tmp_Agents[agent]
                    raise ValueError("Check tmp_agents Sex type %d" % agent)
                else:
                    count_HF += 1
            elif sex_type == 'HM':
                if agent not in self.tmp_HM_agents:
                    print self.tmp_Agents[agent]
                    raise ValueError("Check tmp_agents Sex type %d" % agent)
                else:
                    count_HM += 1
            elif sex_type == 'MSM':
                if agent not in self.tmp_MSM_agents:
                    print self.tmp_Agents[agent]
                    raise ValueError("Check tmp_agents Sex type %d" % agent)
                else:
                    count_MSM += 1
            elif sex_type == 'WSW':
                if agent not in self.tmp_WSW_agents:
                    print self.tmp_Agents[agent]
                    raise ValueError("Check tmp_agents Sex type %d" % agent)
                else:
                    count_WSW += 1
            else:
                raise ValueError("Invalid sex type %s" % str(sex_type))

            # Drug type
            drug_type = agent_dict['Drug Type']
            if drug_type == 'ND':
                if agent not in self.tmp_ND_agents:
                    print self.tmp_Agents[agent]
                    raise ValueError("Check tmp_agents Drug type %d" % agent)
                else:
                    count_ND += 1
            elif drug_type == 'NIDU':
                if agent not in self.tmp_NIDU_agents:
                    print self.tmp_Agents[agent]
                    raise ValueError("Check tmp_agents Drug type %d" % agent)
                else:
                    count_NIDU += 1
            elif drug_type == 'IDU':
                if agent not in self.tmp_IDU_agents:
                    print self.tmp_Agents[agent]
                    raise ValueError("Check tmp_agents Drug type %d" % agent)
                else:
                    count_IDU += 1
            else:
                raise ValueError("Invalid drug type %s" % str(drug_type))

            # HIV
            HIVstatus = agent_dict['HIV']
            if HIVstatus != 0:
                if agent not in self.tmp_HIV_agents:
                    print self.tmp_Agents[agent]
                    raise ValueError("Check tmp_agent HIV %d" % agent)
                else:
                    count_HIV += 1
            # AIDS
            AIDSstatus = agent_dict['AIDS']
            if AIDSstatus != 0:
                if agent not in self.tmp_AIDS_agents:
                    print self.tmp_Agents[agent]
                    raise ValueError("Check agent AIDS %d" % agent)
                else:
                    count_AIDS += 1

        if len(self.tmp_HF_agents) != count_HF:
            raise ValueError("self.tmp_HF agents contains too many agents!")
        if len(self.tmp_HM_agents) != count_HM:
            raise ValueError("self.tmp_HM agents contains too many agents!")
        if len(self.tmp_MSM_agents) != count_MSM:
            raise ValueError("self.tmp_MSM agents contains too many agents!")
        if len(self.tmp_WSW_agents) != count_WSW:
            raise ValueError("self.tmp_WSW agents contains too many agents!")

        if len(self.tmp_NIDU_agents) != count_NIDU:
            raise ValueError("self.tmp_NIDU_agents contains too many agents!")
        if len(self.tmp_ND_agents) != count_ND:
            raise ValueError("self.tmp_ND agents contains too many agents!")
        if len(self.tmp_IDU_agents) != count_IDU:
            mssg = "self.tmp_IDU agents contains too many agents!\
                    \nlen(self.tmp_IDU_agents)=%d\ncount_IDU=%d\n"
            raise ValueError(mssg % (len(self.IDU_agents), count_IDU))

        if len(self.tmp_HIV_agents) != count_HIV:
            raise ValueError("self.tmp_HIV_agents contains too many agents!")
        if len(self.tmp_AIDS_agents) != count_AIDS:
            print "len(self.tmp_AIDS_agents)=%d" % len(self.tmp_AIDS_agents)
            print "count_AIDS=%d" % count_AIDS
            raise ValueError("self.tmp_AIDS agents contains too many agents!")


    def save_AgentPartner_list(self, t):
        """
        :Purpsose:
        Save all agent-partners connections.
        :Input:
        t : int
        Time
        """
        OutFileDir = os.path.expanduser(os.path.join(self.current_dir, 'Results'))
        if not os.path.isdir(OutFileDir):  # create directory if not existing
            os.mkdir(OutFileDir)
        OutFileName = os.path.join(OutFileDir,
                                   'AgentPartnersList_atTime_%s.txt' % str(t))
        if os.path.isfile(OutFileName): os.remove(OutFileName)
        outfile = open(OutFileName, 'w')
        outfile.write('agent\tdrug type\tsex type\tHIV\tAIDS\tHAART\t')
        maxpartners = 0
        for agent in self.Agents:
            numpartners = len(list(self.AdjMat.rows[agent]))
            if numpartners > maxpartners:
                maxpartners = numpartners
        outfile.write('\t'.join(['partner\tp drug type\tp sex type'] *
                                maxpartners))
        outfile.write('\n')
        for agent in sorted(self.Agents.keys()):
            agent_dict = self.Agents[agent]
            outfile.write('%d\t' % agent)
            outfile.write('%s\t' % agent_dict['Drug Type'])
            outfile.write('%s\t' % agent_dict['Sex Type'])
            outfile.write('%d\t' % agent_dict['HIV'])
            outfile.write('%d\t' % agent_dict['AIDS'])
            outfile.write('%d\t' % self.AdherenceAgents[agent])
            for p in sorted(list(self.AdjMat.rows[agent])):
                partner_dict = self.Agents[p]
            outfile.write('%d\t' % int(p))
            outfile.write('%s\t' % partner_dict['Drug Type'])
            outfile.write('%s\t' % partner_dict['Sex Type'])
            outfile.write('\n')


    def _reset_partner_count(self):
        """
        Reset partner count for method assess_interaction_distribution
        """

        # set ND partner count to zero for the next time step
        self.tmp_ND_NumPartners_Count = {}
        self.tmp_NIDU_NumPartners_Count = {}
        self.tmp_IDU_NumPartners_Count = {}
        self.tmp_MSM_NumPartners_Count = {}


    def get_HIV_prevalence_drugs(self):
        """
        get HIV prevalence within all three drug user groups
        """
        count_HIV_IDU = 0
        count_HIV_NIDU = 0
        count_HIV_ND = 0

        for agent in self.Agents:
            HIVstatus = self.get_agent_characteristic(agent, 'HIV')
            if HIVstatus == 1:
                agent_drug_type = self.get_agent_characteristic(agent, 'Drug Type')
            if agent_drug_type == 'IDU':
                count_HIV_IDU += 1
            elif agent_drug_type == 'NIDU':
                count_HIV_NIDU += 1
            elif agent_drug_type == 'ND':
                count_HIV_ND += 1
            elif HIVstatus != 0:
                print HIVstatus
                raise ValueError("HIV status must be either 0 or 1 !")
                # print [count_HIV_IDU, count_HIV_NIDU, count_HIV_ND]
            else:
                raise ValueError("Agent must be either IDU, NIDU or ND !")
        return [count_HIV_IDU, count_HIV_NIDU, count_HIV_ND]


    def get_HIV_prevalence_sex(self):
        """ get HIV prevalence within all four sex groups """
        count_HIV_MSM = 0
        count_HIV_HM = 0
        count_HIV_HF = 0
        count_HIV_WSW = 0

        for agent in self.Agents:
            HIVstatus = self.get_agent_characteristic(agent, 'HIV')
            if HIVstatus == 1:
                agent_sex_type = self.get_agent_characteristic(agent, 'Sex Type')
            if agent_sex_type == 'MSM':
                count_HIV_MSM += 1
            elif agent_sex_type == 'HM':
                count_HIV_HM += 1
            elif agent_sex_type == 'HF':
                count_HIV_HF += 1
            elif agent_sex_type == 'WSW':
                count_HIV_WSW += 1
            elif HIVstatus != 0:
                print HIVstatus
                raise ValueError("HIV status must be either 0 or 1 !")
                # print [count_HIV_IDU, count_HIV_NIDU, count_HIV_ND]
            else:
                raise ValueError("Agent must be either MSM, HM, MF, or WSW !")

        return [count_HIV_MSM, count_HIV_HM, count_HIV_HF, count_HIV_WSW]


    def get_HIV_prevalence_drugs_sex(self):
        """prevalences without and msm only"""
        count_HIV_MIDU = 0
        count_HIV_MNIDU = 0
        count_HIV_MND = 0
        count_HIV_IDUnmsm = 0
        count_HIV_NIDUnmsm = 0
        count_HIV_NDnmsm = 0

        for agent in self.Agents:
            HIVstatus = self.get_agent_characteristic(agent, 'HIV')
            if HIVstatus == 1:
                agent_sex_type = self.get_agent_characteristic(agent, 'Sex Type')
            agent_drug_type = self.get_agent_characteristic(agent, 'Drug Type')
            if agent_drug_type == 'IDU' and agent_sex_type in ['HM', 'HF', 'WSW']:
                count_HIV_IDUnmsm += 1
            elif agent_drug_type == 'IDU' and agent_sex_type == 'MSM':
                count_HIV_MIDU += 1
            elif agent_drug_type == 'NIDU' and agent_sex_type in ['HM', 'HF', 'WSW']:
                count_HIV_NIDUnmsm += 1
            elif agent_drug_type == 'NIDU' and agent_sex_type == 'MSM':
                count_HIV_MNIDU += 1
            elif agent_drug_type == 'ND' and agent_sex_type in ['HM', 'HF', 'WSW']:
                count_HIV_NDnmsm += 1
            elif agent_drug_type == 'ND' and agent_sex_type == 'MSM':
                count_HIV_MND += 1
            elif HIVstatus != 0:
                print HIVstatus
            raise ValueError("HIV status must be either 0 or 1 !")
        return [count_HIV_MIDU, count_HIV_MNIDU, count_HIV_MND, count_HIV_IDUnmsm, count_HIV_NIDUnmsm, count_HIV_NDnmsm]


    def get_HIV_prevalence(self):
        """ get HIV prevalence"""
        HIVcount = 0.0
        for agent in self.Agents.keys():
            HIVstatus = self.get_agent_characteristic(agent, 'HIV')
            if HIVstatus == 1: HIVcount += 1
        return HIVcount


    def return_results(self):
        return self.ResultDict



    def save_result_dict(self):
        OutFileDir = os.path.join(self.current_dir, 'Results')
        if not os.path.isdir(OutFileDir):  # create directory if not existing
            os.mkdir(OutFileDir)
        OutFileName = os.path.join(OutFileDir, 'ResultDictionary.txt')
        if os.path.isfile(OutFileName): os.remove(OutFileName)
        outfile = open(OutFileName, 'w')
        for result_property in sorted(self.ResultDict.keys()):
            outfile.write('%s\t' % result_property)
            for time_t in sorted(self.ResultDict[result_property].keys()):
                outfile.write('%4.5f\t' % float(self.ResultDict[result_property][time_t]))
            outfile.write('\n')


    def save_AdjMat(self, t):
        """
        :Purpose:
        Save Adjacency matrix in sparse format.
        """
        OutFileDir = os.path.expanduser(os.path.join(self.current_dir, 'Results'))
        if not os.path.isdir(OutFileDir):  # create directory if not existing
            os.mkdir(OutFileDir)
        OutFileName = os.path.join(OutFileDir, 'AdjacencyMatrix_atTime_%d.txt' % t)
        if os.path.isfile(OutFileName): os.remove(OutFileName)
        outfile = open(OutFileName, 'w')
        for n, row in enumerate(self.AdjMat.rows):
            outfile.write('%d:\t' % n)
            for partner in row:
                outfile.write('%d,' % partner)
            outfile.write('\n')

    def _writeDNR(self):
        dynnetworkReport = open('Results/dynnetworkReport.txt', 'a')
        for agent in self.Agents:
            sextype = self.Agents[agent]['Sex Type']
            drugtype = self.Agents[agent]['Drug Type']
            HIV = self.Agents[agent]['HIV']
            reportLine = '\t'.join(['0', 'NEWAGENT', repr(agent), sextype, drugtype, repr(HIV)])
            dynnetworkReport.write('\n' + reportLine)
            #open('dynnetworkReport.txt', 'a').write('\n' + reportLine)

        dynnetworkReport.close()
        # initiate HIV status
        for agent in self.Agents:
            if self.Agents[agent]['HIV'] == 1: self._become_HIV(agent, 0)
