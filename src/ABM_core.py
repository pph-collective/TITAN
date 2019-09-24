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

import os
import time
import collections

from scipy.stats import binom
from scipy.stats import poisson
from functools import wraps

try:
    from .HIVABM_Population import PopulationClass, print_population
except ImportError:
    raise ImportError("Can't import PopulationClass")

try:
    from .network_graph_tools import *
except ImportError as e:
    raise ImportError("Can't import network_graph_tools! %s" % str(e))

try:
    from .ABM_partnering import *
except ImportError as e:
    raise ImportError("Can't import ABM_partnering! %s" % str(e))

try:
    from .analysis_output import *
except ImportError as e:
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

        return ret

    return with_profiling


def print_prof_data():
    for fname, data in list(PROF_DATA.items()):
        max_time = max(data[1])
        avg_time = sum(data[1]) / len(data[1])
        print(("Function %s called %d times. " % (fname, data[0])))
        print(
            (
                "\tExecution time max: %.3f, average: %.3f, total %.3f"
                % (max_time, avg_time, sum(data[1]))
            )
        )


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
        returnStr += "Seed: %d\n" % (self.runseed)
        returnStr += "Npop: %d\n" % (params.N_POP)
        returnStr += "Time: %d\n" % (params.TIME_RANGE)
        returnStr += "Mode: %s\n" % (params.model)

        return returnStr

    def __init__(
        self,
        N,
        tmax,
        parameter_dict,
        runseed,
        popseed,
        netseed,
        runtime_diffseed=False,
        model=None,
        network_type=None,
        HIVABM_Agent_set=None,
    ):
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

        if type(tmax) is not int:
            raise ValueError("Number of time steps must be integer")
        else:
            self.tmax = tmax

        if type(runseed) is not int:
            raise ValueError("Random seed must be integer")
        elif runseed == 0:
            self.runseed = random.randint(1, 1000000)
        else:
            self.runseed = runseed

        if type(popseed) is not int:
            raise ValueError("Random seed must be integer")
        elif popseed == 0:
            self.popseed = random.randint(1, 1000000)
        else:
            self.popseed = popseed

        if type(netseed) is not int:
            raise ValueError("Random seed must be integer")
        elif popseed == 0:
            self.netseed = random.randint(1, 1000000)
        else:
            self.netseed = netseed

        self.uniqueSeedID = "r" + str(runseed) + "_p" + str(popseed) + "_n" + str(netseed)

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
            self.Relationships = Agent_set(1, "Relationships")
        else:
            print("\tCreating population Class")
            NetworkClass.__init__(
                self, N=N, network_type=network_type, popSeed=self.popseed, netSeed=self.netseed
            )
            self.All_agentSet.print_subsets()

        self.AdjMat = 0
        self.AdjMats_by_time = 0
        self.newPrEPagents = 0
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
        self.ND_NumPartners = {"ND": [], "NIDU": [], "IDU": [], "MSM": []}  # final counts
        self.NIDU_NumPartners = {"ND": [], "NIDU": [], "IDU": [], "MSM": []}  # final counts
        self.IDU_NumPartners = {"ND": [], "NIDU": [], "IDU": [], "MSM": []}  # final counts
        self.MSM_NumPartners = {"ND": [], "NIDU": [], "IDU": [], "MSM": []}  # final counts

        self.tmp_ND_NumPartners_Count = {}
        self.tmp_NIDU_NumPartners_Count = {}
        self.tmp_IDU_NumPartners_Count = {}
        self.tmp_MSM_NumPartners_Count = {}
        self.tmp_WSW_NumPartners_Count = {}

        self.Acute_agents = []
        self.Transmit_from_agents = []
        self.Transmit_to_agents = []
        self.Transmission_tracker = {"SEX_MSM": {1: 0}, "SEX_NMSM": {1: 0}, "NEEDLE": {1: 0}}
        self.totalDiagnosis = 0
        self.treatmentEnrolled = False

        self.ResultDict = initiate_ResultDict()
        self.newPrEPagents = Agent_set(3, "NewPrEPagents")
        self.newPrEPenrolls = 0
        self.IDUprep = 0
        self.HIVprep = 0
        self.MSMWprep = 0
        # Set seed format. 0: pure random, -1: Stepwise from 1 to nRuns, else: fixed value

        print(("\tRun seed was set to:", runseed))
        self.runRandom = Random(runseed)
        random.seed(self.runseed)
        np.random.seed(self.runseed)
        print(("\tFIRST RANDOM CALL %d" % random.randint(0, 100)))

        print("\tReseting death count")
        self._reset_death_count()  # Number of death

        print("\tCreating network graph")
        self.create_graph_from_agents(self.All_agentSet) #REVIEW redundant with NetworkClass init?

        print("\n === Initialization Protocol Finished ===")

    def run(self, save_adjlist_flag=1, dir_prefix="Results"):
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
            print_stats(
                self,
                self.runseed,
                t,
                self.All_agentSet,
                self.HIV_agentSet,
                self.incarcerated_agentSet,
                self.Trt_PrEP_agentSet,
                self.NewInfections,
                self.NewDiagnosis,
                self.num_Deaths,
                self.ResultDict,
                self.Relationships,
                self.NewHRrolls,
                self.NewIncarRelease,
                self.deathSet,
            )

        def print_components(t):
            name = "componentReport_ALL"
            compReport = open("results/" + name + ".txt", "a")
            components = sorted(nx.connected_component_subgraphs(self.G), key=len, reverse=True)
            compID = 0
            for comp in components:
                totN = nhiv = ntrtmt = ntrthiv = nprep = PrEP_ever_HIV = 0
                for ag in comp.nodes():
                    totN += 1
                    if ag._HIV_bool:
                        nhiv += 1
                        if ag._treatment_bool:
                            ntrthiv += 1
                        if ag._PrEP_ever_bool:
                            PrEP_ever_HIV += 1
                    elif ag._treatment_bool:
                        ntrtmt += 1
                        if ag._PrEP_bool:
                            nprep += 1
                compReport.write(
                    "{rseed}\t{pseed}\t{nseed}\t{t}\t{compID}\t{totalN}\t{Nhiv}\t{Ntrtmt}\t{Nprep}\t{NtrtHIV}\t{NprepHIV}\n".format(
                        rseed=self.runseed,
                        pseed=self.popseed,
                        nseed=self.netseed,
                        t=t,
                        compID=compID,
                        totalN=totN,
                        Nhiv=nhiv,
                        Ntrtmt=ntrtmt,
                        Nprep=nprep,
                        NtrtHIV=ntrthiv,
                        NprepHIV=PrEP_ever_HIV
                    )
                )

                compID += 1
            compReport.close()


        def burnSimulation(burnDuration):
            print(("\n === Burn Initiated for {} timesteps ===".format(burnDuration + 1)))
            for t in range(0, burnDuration + 1):
                self._update_AllAgents(t, burn=True)

                if params.flag_DandR:
                    self._die_and_replace(t)
            print(("\tBurn Cuml Inc:\t{}".format(self.NewInfections.num_members())))
            self.NewInfections.clear_set()
            self.NewDiagnosis.clear_set()
            self.NewHRrolls.clear_set()
            self.NewIncarRelease.clear_set()
            self.newPrEPagents.clear_set()
            self.newPrEPenrolls = 0
            self.IDUprep = 0
            self.HIVprep = 0
            self.MSMWprep = 0

            self._reset_death_count()
            print(" === Simulation Burn Complete ===")

        burnSimulation(params.burnDuration)

        print("\n === Begin Simulation Run ===")
        if params.drawFigures:
            nNodes = self.G.number_of_nodes()
            self.visualize_network(
                coloring=params.drawFigureColor,
                node_size=5000.0 / nNodes,
                curtime=0,
                iterations=10,
                label="Seed" + str(self.runseed),
            )
        if params.calcComponentStats:
            print_components(0)

        self.cumInfT = 0
        self.cumInfW = 0
        self.cumInfB = 0

        def makeAgentZero(numPartners):
            firstHIV = self.runRandom.choice(self.DU_IDU_agentSet._members)
            i = 0
            while i <= numPartners:
                update_partner_assignments(self, 10000.0, self.get_Graph(), agent=firstHIV)
                i += 1
            self._become_HIV(firstHIV, 0)

        print("\t===! Start Main Loop !===")

        # If we are using an agent zero method, create agent zero.
        if params.flag_agentZero:
            makeAgentZero(4)

        if params.drawEdgeList:
            print("Drawing network edge list to file")
            fh = open("results/network/Edgelist_t{}.txt".format(0), "wb")
            self.write_G_edgelist(fh)
            fh.close()
        for t in range(1, self.tmax + 1):
            print(f"\n\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t.: TIME {t}")
            if params.drawFigures and t % params.intermPrintFreq == 0:
                self.visualize_network(
                    coloring=params.drawFigureColor,
                    node_size=5000.0 / nNodes,
                    curtime=t,
                    iterations=10,
                    label="Seed" + str(self.runseed),
                )
            # todo: GET THIS TO THE NEW HIV COUNT

            print(
                (
                    "\tSTARTING HIV count:{}\tTotal Incarcerated:{}\tHR+:{}\tPrEP:{}".format(
                        self.HIV_agentSet.num_members(),
                        self.incarcerated_agentSet.num_members(),
                        self.highrisk_agentsSet.num_members(),
                        self.Trt_PrEP_agentSet.num_members(),
                    )
                )
            )
            print(
                (
                    "Trt:{trt} \t OAT:{oat} \t NAL:{nal}".format(
                        trt=self.treatment_agentSet.num_members(),
                        oat=len(
                            [a for a in self.treatment_agentSet._members if a._OAT_bool == True]
                        ),
                        nal=len(
                            [a for a in self.treatment_agentSet._members if a._naltrex_bool == True]
                        ),
                    )
                )
            )
            self.TimeStep = t

            self._update_AllAgents(t)

            getStats(t)

            self._reset_death_count()

            if params.flag_DandR:
                self._die_and_replace(t)

            print(("Number of relationships: %d" % self.Relationships.num_members()))
            tested = len([tmpA for tmpA in self.HIV_agentSet._members if tmpA._tested])
            self.All_agentSet.print_subsets()

            newInfB = len([tmpA for tmpA in self.NewInfections._members if tmpA._race == "BLACK"])
            newInfW = len([tmpA for tmpA in self.NewInfections._members if tmpA._race == "WHITE"])
            newInfT = len(self.NewInfections._members)
            self.cumInfB += newInfB
            self.cumInfW += newInfW
            self.cumInfT += newInfT

            self.totalDiagnosis += len(self.NewDiagnosis._members)
            if self.totalDiagnosis > params.initTreatment and not self.treatmentEnrolled:
                self._enroll_treatment(t)

            self.NewInfections.clear_set()
            self.NewDiagnosis.clear_set()
            self.NewHRrolls.clear_set()
            self.NewIncarRelease.clear_set()
            self.num_Deaths
            prepReport = open('results/PrEPReport.txt', 'a')
            prepReport.write(
                f"{self.runseed}\t{self.TimeStep}\t{self.newPrEPenrolls}\t{self.IDUprep}\t{self.HIVprep}\t{self.MSMWprep}\n")
            prepReport.close()
            self.newPrEPagents.clear_set()
            self.newPrEPenrolls = 0
            self.IDUprep = 0
            self.HIVprep = 0
            self.MSMWprep = 0
            # If set to draw the edge list, print list at each timestep
            if params.drawEdgeList and t % params.intermPrintFreq == 0:
                print("Drawing network edge list to file")
                fh = open("results/network/Edgelist_t{}.txt".format(t), "wb")
                self.write_G_edgelist(fh)
                fh.close()

            print((t % params.intermPrintFreq))
            if t % params.intermPrintFreq == 0:
                if params.calcNetworkStats:
                    self.write_network_stats(t=t)
                if params.calcComponentStats:
                    print_components(t)

        print_prof_data()


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
        if time == 0:
            i = 0
        elif params.flag_staticN == False:
            update_partner_assignments(self, params.PARTNERTURNOVER, self.get_Graph)
        else:
            pass

        self.Acute_agents = []
        self.Transmit_from_agents = []
        self.Transmit_to_agents = []

        for rel in self.Relationships._members:
            if burn:
                pass
            else:
                self._agents_interact(rel._ID1, rel._ID2, time, rel)
            if params.flag_staticN:
                pass
            else:
                if rel.progress():
                    try:
                        self.get_Graph().remove_edge(rel._ID1, rel._ID2)
                    except:
                        pass
                    self.Relationships.remove_agent(rel)
                    del rel

        if params.flag_HR:
            for tmpA in self.highrisk_agentsSet.iter_agents():
                if tmpA._highrisk_time > 0:
                    tmpA._highrisk_time -= 1
                    if (
                        tmpA._SO == "HM"
                        and params.flag_PrEP
                        and (
                            params.PrEP_target_model == "HR"
                            or params.PrEP_target_model == "IncarHR"
                        )
                    ):
                        for part in tmpA._partners:
                            if not part._HIV_bool:
                                self._initiate_PrEP(part, time)
                else:
                    self.highrisk_agentsSet.remove_agent(tmpA)
                    tmpA._highrisk_bool = False

                    if params.model == "Incar":
                        if tmpA._SO == "HM":
                            tmpA._mean_num_partners -= params.HR_partnerScale
                        elif tmpA._SO == "HF":
                            tmpA._mean_num_partners -= params.HR_partnerScale

        for agent in self.All_agentSet.iter_agents():

            agent_drug_type = agent._DU
            agent_sex_type = agent._SO
            agent_HIV_status = agent._HIV_bool

            agent_incarcerated = agent._incar_bool

            agent._timeAlive += 1
            if params.flag_incar:  # and not burn:
                self._incarcerate(agent, time)

            if agent._MSMW and self.runRandom.random() < params.HIV_MSMW:
                self._become_HIV(agent, 0)

            if agent_drug_type in ["NIDU", "IDU"] and False:
                self._drug_cessation(agent, agent_drug_type)
                self._enter_and_exit_drug_treatment(agent, time)
            if agent_HIV_status:
                if burn:
                    if agent._incar_treatment_time >= 1:
                        agent._incar_treatment_time -= 1

                self._HIVtest(agent, time)
                self._progress_to_AIDS(agent, agent_drug_type)

                if params.flag_ART:
                    self._initiate_HAART(agent, time)
                    agent._HIV_time += 1
            else:
                if params.flag_PrEP:
                    if time >= params.PrEP_startT:
                        if agent._PrEP_bool:
                            self._discont_PrEP(agent, time)
                        elif params.PrEP_target_model == "Clinical":
                            pass
                        elif params.PrEP_target_model == "RandomTrial":
                            pass
                        elif self._PrEP_elligible(agent, time) and not agent._PrEP_bool:
                            self._initiate_PrEP(agent, time)

        if params.flag_PrEP and time >= params.PrEP_startT:
            if params.PrEP_target_model == "Clinical":
                if time > params.PrEP_startT:
                    numPrEP_agents = self.Trt_PrEP_agentSet.num_members()
                    target_PrEP = int(
                        (
                            self.All_agentSet.num_members()
                            - self.All_agentSet._subset["HIV"].num_members()
                        )
                        * params.PrEP_Target
                    )
                    elligiblePool = [
                        ag
                        for ag in self.All_agentSet._subset["SO"]._subset["MSM"]._members
                        if (ag._PrEP_bool == False and ag._HIV_bool == False)
                    ]

                    while numPrEP_agents < target_PrEP:
                        numPrEP_agents = self.Trt_PrEP_agentSet.num_members()
                        target_PrEP = int(
                            (
                                self.All_agentSet.num_members()
                                - self.All_agentSet._subset["HIV"].num_members()
                            )
                            * params.PrEP_Target
                        )
                        self._initiate_PrEP(
                            self._get_clinic_agent(params.PrEP_clinic_cat, elligiblePool), time
                        )
            elif params.PrEP_target_model == "RandomTrial" and time == params.PrEP_startT:
                print("Starting random trial")
                components = sorted(nx.connected_component_subgraphs(self.G), key=len, reverse=True)
                totNods = 0
                for comp in components:
                    totNods += comp.number_of_nodes()
                    if self.runRandom.random() < 0.5:
                        # Component selected as treatment pod!
                        for ag in comp.nodes():
                            if (ag._HIV_bool == False) and (ag._PrEP_bool == False):
                                ag._treatment_bool = True
                                if self.runRandom.random() < params.PrEP_Target:
                                    self._initiate_PrEP(ag, time, force=True)
                print(("Total agents in trial: ", totNods))

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
        partner_HIV_status = agent_2._HIV_bool
        agent_HIV_status = agent_1._HIV_bool
        agent_incar = agent_1._incar_bool
        partner_incar = agent_2._incar_bool
        eligible = False

        # If either agent is incarcerated, skip their interaction
        if agent_incar or partner_incar:
            return

        # Else if neither agent is HIV (shouldn't be possible), skip their interaction to save computation time
        elif not agent_HIV_status and not partner_HIV_status:
            return

        elif agent_HIV_status:  # If agent_1 is HIV
            if partner_HIV_status:  # If agent_1 and agent_2 are both HIV, skip interaction
                return
            else:  # Agent is HIV, partner is succept
                agent = agent_1
                partner = agent_2
                eligible = True
        elif (
            partner_HIV_status
        ):  # If agent_2 is HIV and we have tested both HIV +/-, agent_2 is HIV, agent_1 is succept
            agent = agent_2
            partner = agent_1
            eligible = True

        if eligible:
            partner_drug_type = partner._DU
            agent_drug_type = agent._DU
            partner_sex_type = partner._SO
            agent_sex_type = agent._SO
            partner_HIV_status = partner._HIV_bool
            agent_HIV_status = agent._HIV_bool
            agent_incar = agent._incar_bool
            partner_incar = partner._incar_bool
            if partner_drug_type == "IDU" and agent_drug_type == "IDU":
                # Injection is possible
                # If agent is on post incar HR treatment to prevent IDU behavior, pass IUD infections
                if agent._incar_treatment_time > 0 and params.inc_treat_IDU_beh:
                    pass

                elif self._sex_possible(agent_sex_type, partner_sex_type):
                    # Sex is possible
                    rv = self.runRandom.random()
                    if rv < 0.25:  # Needle only (60%)
                        self._needle_transmission(agent, partner, time)
                    else:  # Both sex and needle (20%)
                        self._needle_transmission(agent, partner, time)
                        self._sex_transmission(agent, partner, time, rel)
                else:
                    # Sex not possible, needle only
                    self._needle_transmission(agent, partner, time)

            elif partner_drug_type in ["NIDU", "NDU"] or agent_drug_type in ["NIDU", "NDU"]:
                if self._sex_possible(agent_sex_type, partner_sex_type):
                    self._sex_transmission(agent, partner, time, rel)
                else:
                    return
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

        partner_drug_type = self.get_agent_characteristic(partner, "Drug Type")
        agent_drug_type = self.get_agent_characteristic(agent, "Drug Type")
        Flag_Partner_IDU_NIDU_Transition = 0
        Flag_Agent_IDU_NIDU_Transition = 0

        # NIDU -> IDU
        if agent_drug_type == "NIDU" and partner_drug_type == "IDU":
            if self.runRandom.random() < 0.00875 / 12:
                self.tmp_Agents[agent].update({"Drug Type": "IDU"})  # agent becomes IDU
            # Sex type lists
            if agent in self.tmp_NIDU_agents:
                self.tmp_NIDU_agents.remove(agent)
            if agent in self.tmp_ND_agents:  # agent might have transitioned into ND before
                self.tmp_ND_agents.remove(agent)
            if agent not in self.tmp_IDU_agents:
                self.tmp_IDU_agents.append(agent)

        elif partner_drug_type == "NIDU" and agent_drug_type == "IDU":
            if self.runRandom.random() < 0.0175 / 12:
                self.tmp_Agents[partner].update({"Drug Type": "IDU"})  # partner becomes IDU
            # Sex type lists
            if partner in self.tmp_NIDU_agents:
                self.tmp_NIDU_agents.remove(partner)
            if partner in self.tmp_ND_agents:  # agent might have transitioned into ND before
                self.tmp_ND_agents.remove(partner)
            if partner not in self.tmp_IDU_agents:
                self.tmp_IDU_agents.append(partner)

        ## ND -> IDU
        elif agent_drug_type == "ND" and partner_drug_type == "IDU":
            if self.runRandom.random() < 0.001:
                self.tmp_Agents[agent].update({"Drug Type": "IDU"})  # agent becomes NIDU
            # Sex type lists
            if agent in self.tmp_ND_agents:  # agent might have transitioned already
                self.tmp_ND_agents.remove(agent)
            if agent not in self.tmp_IDU_agents:
                self.tmp_IDU_agents.append(agent)

        # ND -> IDU
        elif partner_drug_type == "ND" and agent_drug_type == "IDU":
            if self.runRandom.random() < 0.001:
                self.tmp_Agents[agent].update({"Drug Type": "IDU"})  # agent becomes NIDU
            # Sex type lists
            if agent in self.tmp_ND_agents:  # agent might have transitioned already
                self.tmp_ND_agents.remove(agent)
            if agent not in self.tmp_IDU_agents:
                self.tmp_IDU_agents.append(agent)

        # ND -> NIDU
        elif agent_drug_type == "ND" and partner_drug_type == "NIDU":
            if self.runRandom.random() < 0.005:
                self.tmp_Agents[agent].update({"Drug Type": "NIDU"})  # agent becomes NIDU
            # Sex type lists
            if agent in self.tmp_ND_agents:  # agent might have transitioned already
                self.tmp_ND_agents.remove(agent)
            if agent not in self.tmp_NIDU_agents:
                self.tmp_NIDU_agents.append(agent)
        # ND -> NIDU
        elif partner_drug_type == "ND" and agent_drug_type == "NIDU":
            if self.runRandom.random() < 0.005:
                self.tmp_Agents[partner].update({"Drug Type": "NIDU"})  # partner becomes NIDU
            # Sex type lists
            if partner in self.tmp_ND_agents:  # agent might have transitioned already
                self.tmp_ND_agents.remove(partner)
            if partner not in self.tmp_NIDU_agents:
                self.tmp_NIDU_agents.append(partner)

        # NIDU -> ND from agent's perspective
        elif agent_drug_type == "NIDU" and partner_drug_type == "ND":
            if self.runRandom.random() < 0.001:
                self.tmp_Agents[agent].update({"Drug Type": "ND"})
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
        elif partner_drug_type == "NIDU" and agent_drug_type == "ND":
            if self.runRandom.random() < 0.001:
                self.tmp_Agents[partner].update({"Drug Type": "ND"})
            # Sex type lists
            if partner in self.tmp_NIDU_agents:  # partner might have transitioned already
                self.tmp_NIDU_agents.remove(partner)
            if (
                partner in self.tmp_IDU_agents
            ):  # partner might have previously transitioned into IDU
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
        acuteTimePeriod = 2
        hiv_t = agent._HIV_time

        if hiv_t <= acuteTimePeriod and hiv_t > 0:
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

        sex_type = agent._SO
        race_type = agent._race
        ageBin = agent._ageBin
        tested = agent._tested
        onHAART = agent._HAART_bool

        agentAdherence = agent._HAART_adh
        "Logic for if needle or sex type interaction"
        if interaction == "NEEDLE":
            p = params.TransmissionProbabilities["NEEDLE"][str(agentAdherence)]

        elif interaction == "SEX":
            p = params.TransmissionProbabilities["SEX"][sex_type][str(agentAdherence)]

        isAcute = self.get_acute_status(agent, 0)

        # Scaling parameter for acute HIV infections
        if isAcute:
            p = p * params.cal_AcuteScaling

        # Scaling parameter for positively identified HIV agents
        if tested:
            p = p * (1 - params.cal_RR_Dx)

        # Tuning parameter for ART efficiency
        if onHAART:  # self.AdherenceAgents[agent] > 0:
            p = p * params.cal_RR_HAART

        # Racial calibration parameter to attain proper race incidence disparity
        if race_type == "BLACK":
            p = p * params.cal_raceXmission

        # Scaling parameter for per act transmission.
        p = p * params.cal_pXmissionScaling

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
        # Param to scale number of partners

        # both must be IDU
        partner_drug_type = partner._DU
        agent_drug_type = agent._DU
        agent_race = agent._race
        agent_sex_type = agent._SO
        Race_Agent = agent._race
        Type_agent = agent._SO

        if not (partner_drug_type == "IDU" and agent_drug_type == "IDU"):
            raise ValueError(
                "To share a needle both agents must be IDU!%s %s"
                % (str(agent_drug_type), str(partner_drug_type))
            )
        NumberP = len(self.ExistingLinksCollapsedList)
        # Do they share a needle?
        SEPstat = agent._SNE_bool

        isAcute = self.get_acute_status(agent, time)
        HIV_agent = agent._HIV_bool
        HIV_partner = partner._HIV_bool
        MEAN_N_ACTS = (
            params.DemographicParams[Race_Agent][Type_agent]["NUMSexActs"]
            * params.cal_NeedleActScaling
        )
        share_acts = poisson.rvs(MEAN_N_ACTS, size=1)[0]

        if SEPstat:
            p_UnsafeNeedleShare = 0.02  # no needle sharing
        else:  # they do share a needle
            # HIV+ ?
            if share_acts < 1:
                share_acts = 1

            p_UnsafeNeedleShare = (
                params.DemographicParams[agent_race][agent_sex_type]["NEEDLESH"]
                * params.safeNeedleExchangePrev
            )

        for n in range(share_acts):
            if self.runRandom.random() > p_UnsafeNeedleShare:
                share_acts -= 1

        if HIV_agent == 1 and HIV_partner == 0 and share_acts >= 1.0:
            p = self.get_transmission_probability(agent, "NEEDLE")
            p_transmission = binom.pmf(1.0, share_acts, p)

            p_total_transmission = 0
            if share_acts == 1:
                p_total_transmission = p
            else:
                p_total_transmission = 1.0 - binom.pmf(0, share_acts, p)

            if self.runRandom.random() < p_total_transmission:
                # if agent HIV+ partner becomes HIV+
                self._become_HIV(partner, time)
                self.Transmit_from_agents += [agent]
                self.Transmit_to_agents += [partner]
                if isAcute:
                    self.Acute_agents += [agent]
                else:
                    pass

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
        # Double check: Sex possible?
        Type_agent = agent._SO
        Type_partner = partner._SO
        if not self._sex_possible(Type_agent, Type_partner):
            raise ValueError("Sex must be possible! %s %s" % (str(Type_agent), str(Type_partner)))

        # HIV status of agent and partner
        # Everything from here is only run if one of them is HIV+

        HIVstatus_Agent = agent._HIV_bool
        HIVstatus_Partner = partner._HIV_bool
        AIDSstatus_Agent = agent._AIDS_bool
        AIDSstatus_Partner = partner._AIDS_bool
        Race_Agent = agent._race
        isAcute = self.get_acute_status(agent, time)

        if HIVstatus_Partner:
            pass
        if HIVstatus_Agent and HIVstatus_Partner:
            return
        elif HIVstatus_Agent == 1 or HIVstatus_Partner == 1:
            # Sex between men?
            if Type_agent == "MSM" and Type_partner == "MSM":
                SexBetweenMen = 1
            else:
                SexBetweenMen = 0
            # Define probabilities for unsafe sex

            # unprotected sex probabilities for primary partnerships
            p_UnsafeSafeSex1 = params.DemographicParams[Race_Agent][Type_agent]["UNSAFESEX"]
            MSexActs = self._get_number_of_sexActs(agent) * params.cal_SexualActScaling
            T_sex_acts1 = int(poisson.rvs(MSexActs, size=1))

            num_int = rel._total_sex_acts
            # Get condom usage
            if params.condomUseType == "Race":
                p_UnsafeSafeSex1 = params.DemographicParams[Race_Agent][Type_agent]["UNSAFESEX"]
            else:
                if num_int < 10:
                    if num_int == 0:
                        p_UnsafeSafeSex1 = 0.443
                    elif num_int == 1:
                        p_UnsafeSafeSex1 = 0.481
                    else:
                        p_UnsafeSafeSex1 = 0.514
                else:  # More than 10 acts
                    p_UnsafeSafeSex1 = 0.759

            # Reduction of risk acts between partners for condom usage
            U_sex_acts1 = T_sex_acts1
            for n in range(U_sex_acts1):
                if self.runRandom.random() < p_UnsafeSafeSex1:
                    U_sex_acts1 -= 1

            U_sex_acts2 = U_sex_acts1

            if U_sex_acts2 >= 1:
                # if agent HIV+
                rel._total_sex_acts += U_sex_acts2
                if HIVstatus_Agent == 1 or HIVstatus_Partner == 1:
                    ppAct = self.get_transmission_probability(agent, "SEX")

                    # Reduction of transmissibility for acts between partners for PrEP adherence
                    if agent._PrEP_bool or partner._PrEP_bool:
                        if agent._PrEPresistance or partner._PrEPresistance:
                            pass

                        elif params.PrEP_type == "Oral":
                            if agent._PrEP_adh == 1 or partner._PrEP_adh == 1:
                                ppAct = ppAct * (1.0 - params.PrEP_AdhEffic)  # 0.04
                            else:
                                ppAct = ppAct * (1.0 - params.PrEP_NonAdhEffic)  # 0.24

                        elif params.PrEP_type == "Inj":
                            ppActReduction = -1.0 * np.exp(-5.528636721 * partner._PrEP_load) + 1
                            if agent._PrEP_adh == 1 or partner._PrEP_adh == 1:
                                ppAct = ppAct * (1.0 - ppActReduction)  # 0.04

                    p_total_transmission = 0
                    if U_sex_acts2 == 1:
                        p_total_transmission = ppAct
                    else:
                        p_total_transmission = 1.0 - binom.pmf(0, U_sex_acts1, ppAct)

                    if self.runRandom.random() < p_total_transmission:

                        # if agent HIV+ partner becomes HIV+
                        self.Transmit_from_agents += [agent]
                        self.Transmit_to_agents += [partner]
                        self._become_HIV(partner, time)
                        if isAcute:
                            self.Acute_agents += [agent]
                        else:
                            pass
            else:
                return


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
        i = 0

        while True:
            i += 1
            pMatch += params.sexualFrequency[i]["p_value"]
            if rv <= pMatch:
                minSA = params.sexualFrequency[i]["min"]
                maxSA = params.sexualFrequency[i]["max"]
                return self.runRandom.randrange(minSA, maxSA, 1)
            if i == 5:
                break


    def _become_HIV(self, agent, time):
        """
        :Purpose:
            agent becomes HIV agent. Update all appropriate list and
            dictionaries.

        :Input:
            agent : int
            time : int

        """
        if agent._HIV_bool:
            pass
        else:
            agent._HIV_bool = True
            agent._HIV_time = 1
            self.NewInfections.add_agent(agent)
            self.HIV_agentSet.add_agent(agent)
            if agent._PrEP_time > 0:
                if self.runRandom.random() < params.PrEP_resist:
                    agent._PrEPresistance = 1

        if agent._PrEP_bool:
            self._discont_PrEP(agent, time, force=True)


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
            raise ValueError("Invalid partner_sex_type! %s" % str(partner_sex_type))

        # Sex possible
        if agent_sex_type == "HM" and partner_sex_type in ["HF", "WSW", "MTF"]:
            SexPossible = True
        elif agent_sex_type == "MSM" and partner_sex_type in ["MSM", "WSW", "HF", "MTF"]:
            SexPossible = True
        elif agent_sex_type == "WSW" and partner_sex_type in ["MSM", "WSW", "HM"]:
            SexPossible = True
        elif agent_sex_type == "HF" and partner_sex_type in ["HM", "MSM"]:
            SexPossible = True
        elif agent_sex_type == "MTF" and partner_sex_type in ["HM", "MSM"]:
            SexPossible = True
        else:
            SexPossible = False

        if agent_sex_type == "HM" and partner_sex_type == "HM" and SexPossible: #REVIEW move to tests
            raise ValueError("Check _sex_possible method!")

        return SexPossible

    def _drug_cessation(self, agent, agent_drug_type):
        """
        :Purpose:
            Account for drug cessation of IDU to NIDU and NIDU to ND.

        :Input:
            agent : int
            agent_drug_type : str ("IDU", "NIDU")

        """
        if agent_drug_type == "IDU":
            if self.runRandom.random() < 0.017 / 12 / 2:
                self.tmp_Agents[agent].update({"Drug Type": "NIDU"})
                if agent in self.tmp_IDU_agents:  # agent might have transitioned already
                    self.tmp_IDU_agents.remove(agent)
                if agent not in self.tmp_NIDU_agents:
                    self.tmp_NIDU_agents.append(agent)
        elif agent_drug_type == "NIDU":
            if self.runRandom.random() < 0.017 / 12:
                self.tmp_Agents[agent].update({"Drug Type": "ND"})
                if agent in self.tmp_NIDU_agents:  # agent might have transitioned already
                    self.tmp_NIDU_agents.remove(agent)
                if agent not in self.tmp_ND_agents:
                    self.tmp_ND_agents.append(agent)
                if agent in self.DrugTreatmentAgents_current:
                    self._exit_drug_treatment(agent)
        else:
            raise ValueError("Drug cessation only valid for IDU and NIDU!")

    def _enroll_treatment(self, time):
        """
        :Purpose:
            Account for drug cessation of IDU to NIDU and NIDU to ND.

        :Input:
            time : int

        """
        print(("\n\n!!!!Engaginge treatment process: %d" % time))
        self.treatmentEnrolled = True
        for agent in self.All_agentSet.iter_agents():
            if self.runRandom.random() < params.treatmentCov and agent._DU == "IDU":
                agent._SNE_bool = True


    def _becomeHighRisk(self, agent, HRtype=None, duration=None):

        if agent not in self.highrisk_agentsSet._members:
            self.highrisk_agentsSet.add_agent(agent)
        if not agent._everhighrisk_bool:
            self.NewHRrolls.add_agent(agent)

        agent._highrisk_bool = True
        agent._everhighrisk_bool = True
        agent._highrisk_type = HRtype

        if duration != None:
            agent._highrisk_time = duration
        else:
            agent._highrisk_time = params.HR_M_dur

    def _incarcerate(self, agent, time):
        """
        :Purpose:
            Account for drug cessation of IDU to NIDU and NIDU to ND.

        :Input:
            agent : int

        """
        drug_type = agent._DU
        sex_type = agent._SO
        race_type = agent._race
        hiv_bool = agent._HIV_bool
        tested = agent._tested
        incar_t = agent._incar_time
        incar_bool = agent._incar_bool
        haart_bool = agent._HAART_bool

        if incar_bool:
            agent._incar_time -= 1

            # get out if t=0
            if incar_t == 1:  # FREE AGENT
                self.incarcerated_agentSet.remove_agent(agent)
                self.NewIncarRelease.add_agent(agent)
                agent._incar_bool = False
                agent._ever_incar_bool = True
                if params.model == "Incar":
                    if (
                        not agent._highrisk_bool and params.flag_HR
                    ):  # If behavioral treatment on and agent HIV, ignore HR period.
                        if (
                            params.inc_treat_HRsex_beh
                            and hiv_bool
                            and (time >= params.inc_treatment_startdate)
                        ):
                            pass
                        else:  # Else, become high risk
                            self.highrisk_agentsSet.add_agent(agent)
                            if not agent._everhighrisk_bool:
                                self.NewHRrolls.add_agent(agent)

                            agent._mean_num_partners = (
                                agent._mean_num_partners + params.HR_partnerScale
                            )
                            agent._highrisk_bool = True
                            agent._everhighrisk_bool = True
                            agent._highrisk_time = params.HR_M_dur

                    if (
                        params.inc_treat_RIC
                        or params.inc_treat_HRsex_beh
                        or params.inc_treat_IDU_beh
                    ) and (time >= params.inc_treatment_startdate):
                        agent._incar_treatment_time = params.inc_treatment_dur

                    if hiv_bool:
                        if haart_bool:
                            if (
                                self.runRandom.random() > params.inc_ARTdisc
                            ):  # 12% remain surpressed
                                pass

                            else:
                                agent._HAART_bool = False
                                agent._HAART_adh = 0
                                self.Trt_ART_agentSet.remove_agent(agent)

                            ### END FORCE ####
                elif params.model == "Overdose":
                    self._becomeHighRisk(agent, HRtype="postIncar", duration=26)
                    if self.runRandom.random() < params.p_enroll_OAT_post_release:
                        self._enter_drug_treatment(agent, trtType="OAT")
                        self._DOC_OAT_bool = True
                    elif self.runRandom.random() < params.p_enroll_Nal_post_release:
                        self._enter_drug_treatment(agent, trtType="NAL")
                        self._DOC_NAL_bool = True

        elif (
            self.runRandom.random()
            < params.DemographicParams[race_type][sex_type]["INCAR"]
            * (1 + (hiv_bool * 4))
            * params.cal_IncarP
        ):
            toss = 2
            # M 61% for six months or less, 13% for six months to one year, 16% for 1-3 years, 4% for 3-5 years, 4% for 5-10 years
            # F 22% for 6 months or less, 17% for 6mo-1 year, 24% for 1-3 years, 14% for 3-5 years 13% for 5-10years, 10% for 10+ years
            if agent._SO == "HF":
                jailDuration = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}, 6: {}}
                jailDuration[1] = {"p_value": (0.40), "min": 1, "max": 2}
                jailDuration[2] = {"p_value": (0.475), "min": 1, "max": 13}
                jailDuration[3] = {"p_value": (0.065), "min": 13, "max": 26}
                jailDuration[4] = {"p_value": (0.045), "min": 26, "max": 78}
                jailDuration[5] = {"p_value": (0.01), "min": 78, "max": 130}
                jailDuration[6] = {"p_value": (0.01), "min": 130, "max": 260}

            elif agent._SO == "HM":
                jailDuration = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}, 6: {}}
                jailDuration[1] = {"p_value": (0.43), "min": 1, "max": 2}
                jailDuration[2] = {"p_value": (0.50), "min": 1, "max": 13}
                jailDuration[3] = {"p_value": (0.02), "min": 13, "max": 26}
                jailDuration[4] = {"p_value": (0.02), "min": 26, "max": 78}
                jailDuration[5] = {"p_value": (0.03), "min": 78, "max": 130}
                jailDuration[6] = {"p_value": (0.01), "min": 130, "max": 260}
            if toss == 1:  # JAIL
                timestay = self.runRandom.randint(
                    params.inc_JailMin, params.inc_JailMax
                )
                if hiv_bool and not tested:
                    if self.runRandom.random() < params.inc_JailTestProb:
                        agent._tested = True

            # M 13% for 6 months or less, 8% for 6 mo-1year, 20% for 1-3 years, 11% for 3-5 years, 16% for 5-10 years, 30% for 10+ years
            # F 74% for 6 months or less, 12% for 6 months to one year, 10% for one to three years, 4% for over three years
            else:  # PRISON
                durationBin = current_p_value = 0
                p = self.runRandom.random()
                while p > current_p_value:
                    durationBin += 1
                    current_p_value += jailDuration[durationBin]["p_value"]
                timestay = self.runRandom.randint(
                    jailDuration[durationBin]["min"], jailDuration[durationBin]["max"]
                )

                if hiv_bool:
                    if not tested:
                        if self.runRandom.random() < params.inc_PrisTestProb:
                            agent._tested = True
                    else:  # Then tested and HIV, check to enroll in ART
                        if self.runRandom.random() < params.inc_ARTenroll:
                            tmp_rnd = self.runRandom.random()
                            HAART_ADH = params.inc_ARTadh
                            if tmp_rnd < HAART_ADH:
                                adherence = 5
                            else:
                                adherence = self.runRandom.randint(1, 4)

                            # Add agent to HAART class set, update agent params
                            agent._HAART_bool = True
                            agent._HAART_adh = adherence
                            agent._HAART_time = time
                            self.Trt_ART_agentSet.add_agent(agent)

            agent._incar_bool = True
            agent._incar_time = timestay
            self.incarcerated_agentSet.add_agent(agent)
            self.totalIncarcerated += 1

            # PUT PARTNERS IN HIGH RISK
            for tmpA in agent._partners:
                if tmpA._highrisk_bool == True:
                    pass
                else:
                    if self.runRandom.random() < params.HR_proportion:
                        if not tmpA._highrisk_bool:
                            self.highrisk_agentsSet.add_agent(tmpA)
                            if not tmpA._everhighrisk_bool:
                                self.NewHRrolls.add_agent(tmpA)
                            tmpA._mean_num_partners += (
                                params.HR_partnerScale
                            )  # 32.5 #2 + 3.25 from incar HR
                            tmpA._highrisk_bool = True
                            tmpA._everhighrisk_bool = True
                            tmpA._highrisk_time = params.HR_F_dur
                if params.flag_PrEP and (
                    params.PrEP_target_model == "Incar" or params.PrEP_target_model == "IncarHR"
                ):
                    # Atempt to put partner on prep if less than probability
                    if not tmpA._HIV_bool:
                        self._initiate_PrEP(tmpA, time)

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

        drug_type = agent._DU
        sex_type = agent._SO
        race_type = agent._race
        hiv_Status = agent._HIV_bool
        tested = agent._tested
        if not tested:
            test_prob = params.DemographicParams[race_type][sex_type]["HIVTEST"]

            # Rescale based on calibration param
            test_prob = test_prob * params.cal_TestFreq

            # If roll less than test probablity
            if self.runRandom.random() < test_prob:  ###WAS / 10
                # Become tested, add to tested agent set
                agent._tested = True
                self.NewDiagnosis.add_agent(agent)
                self.Trt_Tstd_agentSet.add_agent(agent)
                # If treatment co-enrollment enabled and coverage greater than 0
                if self.treatmentEnrolled and params.treatmentCov > 0:
                    # For each partner, attempt to test for HIV
                    for ptnr in agent._partners:
                        if ptnr._HIV_bool and not ptnr._tested:
                            if self.runRandom.random() < 0.87:
                                ptnr._tested = True
                                self.NewDiagnosis.add_agent(ptnr)

    #REVIEW not used anywhere
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
        drug_type = self.get_agent_characteristic(agent, "Drug Type")
        SEPstat = False
        if agent in self.SEPAgents:
            if time == self.SEPAgents[agent]:
                SEPstat = True
        if drug_type == "IDU":
            if SEPstat:
                if self.runRandom.random() < self.VCT_NSP:  # !!!!!!!!!!!!!!!!!!!!
                    self.VCTAgents.update({agent: time})
            else:
                if self.runRandom.random() < self.VCT_NoNSP_IDU:  # !!!!!!!!!!!!!!!!!!!
                    self.VCTAgents.update({agent: time})
        if drug_type == "NIDU":
            if self.runRandom.random() < self.VCT_NoNSP_NIDU:  # !!!!!!!!!!!!!!!!!!
                self.VCTAgents.update({agent: time})
        elif agent in self.MSM_agents and self.runRandom.random() < self.VCT_NoNSP_MSM:  # !
            self.VCTAgents.update({agent: time})
        else:
            if self.runRandom.random() < self.VCT_NoNSP_EE:  # !!!!!!!!!!!!!!!!!!!!!!!!!
                self.VCTAgents.update({agent: time})

    #REVIEW not used anywhere
    def _SEP(self, agent, time):
        """
        :Purpose:
            Account for SEP (Syringe Exchange Program) for IDU agents. \n

        :Input:
            agent : Agent
            time : int

        :Output:
            SEPstat : bool
        """
        if self.get_agent_characteristic(agent, "Drug Type") != "IDU":
            raise ValueError(
                "_SEP only valid for IDU agents! agent: %d %s"
                % (agent, str(self.get_agent_characteristic(agent, "Drug Type")))
            )
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

    off = 0


    def _exit_drug_treatment(self, agent):
        """
        Agent exits drug treament.
        """

        if agent._OAT_bool:
            agent._OAT_bool = False
            self._becomeHighRisk(agent, HRtype="postTrtOAT", duration=1)
            agent._off = True
        if agent._naltrex_bool:
            agent._naltrex_bool = False
            self._becomeHighRisk(agent, HRtype="postTrtNal", duration=1)
            agent._off = True

        agent._treatment_bool = False
        self._naltrex_bool = False
        self._OAT_bool = False
        agent._treatment_time = 0
        self.treatment_agentSet.remove_agent(agent)


    def _enter_drug_treatment(self, agent, trtType=None):
        agent._treatment_bool = True

        if trtType == "OAT":
            agent._OAT_bool = True
        elif trtType == "NAL":
            agent._naltrex_bool = True
        elif self.runRandom.random() < params.MATasOAT:
            agent._OAT_bool = True
        else:
            agent._naltrex_bool = True

        try:
            self.treatment_agentSet.add_agent(agent)
        except:
            print(
                (
                    "agent %s is already a member of agent set %s"
                    % (agent.get_ID(), targetSet.get_ID())
                )
            )

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
        # SEB edits here
        agent_race = agent._race
        agent_so = agent._SO
        MATProb = params.DemographicParams[agent_race][agent_so]["MATProbScalar"]
        discMATProb = params.DemographicParams[agent_race][agent_so]["MAT_disc_prob"]
        if agent._treatment_bool:
            if self.runRandom.random() < discMATProb:
                self._exit_drug_treatment(agent)
        elif self.runRandom.random() < MATProb:
            self._enter_drug_treatment(agent)
        else:
            pass


    def _initiate_HAART(self, agent, time):
        """
        :Purpose:
            Account for HIV treatment through highly active antiretroviral therapy (HAART).
            HAART was implemented in 1996, hence, there is treatment only after 1996.
            HIV treatment assumes that the agent knows their HIV+ status.

        :Input:
            agent : Agent
            time : int

        :Output:
            none
        """

        HAART_coverage = 0.30

        # Check valid input
        if not agent._HIV_bool:
            print(("HIV_agents: ", sorted(self.HIV_agents)))
            print(("tmp_HIV_agents: ", sorted(self.tmp_HIV_agents)))
            print(("Agent[agent]", self.Agents[agent]))
            try:
                print(("tmp_Agent[agent]", self.tmp_Agents[agent]))
            except KeyError:
                pass
            raise ValueError("HAART only valid for HIV agents!agent:%s" % str(agent))

        agent_drug_type = agent._DU
        agent_haart = agent._HAART_bool
        agent_HIV = agent._HIV_bool
        agent_Test_bool = agent._tested
        agent_race = agent._race
        agent_so = agent._SO

        # Set HAART agents adherence at t=0 for all instanced HAART
        if time == 0 and agent_haart:
            agent_haart_adh = agent._HAART_adh
            if agent_haart_adh == 0:
                tmp_rnd = self.runRandom.random()
                HAART_ADH = params.DemographicParams[agent_race][agent_so]["HAARTadh"]
                if tmp_rnd < HAART_ADH:
                    adherence = 5
                else:
                    adherence = self.runRandom.randint(1, 4)

                # add to agent haart set
                agent._HAART_adh = adherence
                agent._HAART_time = time

        # Determine probability of HIV treatment
        if time >= 0 and agent_Test_bool:
            # Go on HAART
            if not agent_haart and agent._HAART_time == 0:
                if self.runRandom.random() < params.DemographicParams[agent_race][agent_so]["HAARTprev"] * params.cal_ART_cov:

                    tmp_rnd = self.runRandom.random()
                    HAART_ADH = params.DemographicParams[agent_race][agent_so]["HAARTadh"]
                    if tmp_rnd < HAART_ADH:
                        adherence = 5
                    else:
                        adherence = self.runRandom.randint(1, 4)

                    # Add agent to HAART class set, update agent params
                    agent._HAART_bool = True
                    agent._HAART_adh = adherence
                    agent._HAART_time = time
                    self.Trt_ART_agentSet.add_agent(agent)
            elif (
                agent_haart
                and self.runRandom.random()
                < params.DemographicParams[agent_race][agent_so]["HAARTdisc"]
            ):
                if agent._incar_treatment_time > 0 and params.inc_treat_RIC:
                    pass
                else:
                    agent._HAART_bool = False
                    agent._HAART_adh = 0
                    agent._HAART_time = 0
                    self.Trt_ART_agentSet.remove_agent(agent)

    def _PrEP_elligible(self, agent, time):
        elligble = False
        if params.PrEP_target_model == "Allcomers":
            elligble = True
        elif params.PrEP_target_model == "CDCwomen":
            if agent._SO == 'HF':
                for ptn in set(agent._relationships):
                    if ptn._ID1 == agent:
                        partner = ptn._ID2
                    else:
                        partner = ptn._ID1
                    if ptn._duration > 1:
                        if partner._DU == 'IDU':
                            elligble = True
                            agent._PrEP_reason.append('IDU')
                        if partner._tested:
                            elligble = True
                            agent._PrEP_reason.append('HIV test')
                        if partner._MSMW:
                            elligible = True
                            agent._PrEP_reason.append('MSMW')
        elif params.PrEP_target_model == "CDCmsm":
            if agent._SO == "MSM":
                for ptn in agent._relationships:
                    if ptn._ID1 == agent:
                        partner = ptn._ID2
                    else:
                        partner = ptn._ID1
                    if ptn._duration > 1:
                        if (partner._tested or agent._mean_num_partners > 1):
                            elligble = True
                            break
        elif params.PrEP_target_model == "HighPN5":
            if agent._mean_num_partners >= 5:
                elligble = True
        elif params.PrEP_target_model == "HighPN10":
            if agent._mean_num_partners >= 10:
                elligble = True
        elif params.PrEP_target_model == "SRIns":
            if agent._sexualRole == "Insertive":
                elligble = True
        elif params.PrEP_target_model == "MSM":
            if agent._SO == ("MSM" or "MTF"):
                elligble = True
        elif params.PrEP_target_model == "RandomTrial":
            # If using random trial
            if time == 0:
                # if in init timestep 0, use agent set elligiblity
                elligible = agent._PrEP_elligible
            if time > 0:
                # else, false to not continue enrollment past random trial start
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
        # N(t) = N0 (0.5)^(t/t_half)
        agent._PrEP_lastDose += 1
        if agent._PrEP_lastDose > 12:
            agent._PrEP_load = 0.0
        else:
            agent._PrEP_load = params.PrEP_peakLoad * (
                (0.5) ** (agent._PrEP_lastDose / (params.PrEP_halflife / 30))
            )


    def _discont_PrEP(self, agent, time, force=False):
        # If force flag set, auto kick off prep.
        if force == True:
            self.Trt_PrEP_agentSet.remove_agent(agent)
            agent._PrEP_bool = False
            agent._PrEP_reason = []
        # else if agent is no longer enrolled on PrEP, increase time since last dose
        elif agent._PrEP_time > 0:
            if agent._PrEP_time == 1:
                agent._PrEP_bool = False
                agent._PrEP_reason = []
                agent._PrEP_time -= 1
            else:
                agent._PrEP_time -= 1

        # else if agent is on PrEP, see if they should discontinue
        elif agent._PrEP_bool and agent._PrEP_time == 0:
            if (
                self.runRandom.random()
                < params.DemographicParams[agent._race][agent._SO]["PrEPdisc"]
            ):
                agent._PrEP_time = params.PrEP_falloutT
                self.Trt_PrEP_agentSet.remove_agent(agent)

                if params.PrEP_type == "Oral":
                    agent._PrEP_bool = False
                    agent._PrEP_reason = []
            else:  # if not discontinue, see if its time for a new shot.
                if agent._PrEP_lastDose > 2:
                    agent._PrEP_lastDose = -1

        if params.PrEP_type == "Inj":
            self._calc_PrEP_load(agent)

    def _initiate_PrEP(self, agent, time, force=False):
        """
        :Purpose:
            Place agents onto PrEP treatment.
            PrEP treatment assumes that the agent knows his HIV+ status is negative.

        :Input:
            agent : Agent
            time : int
            force : default is `False`

        :Output:
            none
        """

        def _enrollPrEP(self, agent):
            agent._PrEP_bool = True
            agent._PrEP_time = 0
            self.Trt_PrEP_agentSet.add_agent(agent)
            self.newPrEPagents.add_agent(agent)
            self.newPrEPenrolls += 1
            if params.PrEP_target_model == "CDCwomen":
                if 'IDU' in agent._PrEP_reason:
                    self.IDUprep += 1
                if 'HIV test' in agent._PrEP_reason:
                    self.HIVprep += 1
                if 'MSMW' in agent._PrEP_reason:
                    self.MSMWprep += 1
            tmp_rnd = self.runRandom.random()
            if params.PrEP_Adherence == 'AtlantaMSM':
                if tmp_rnd < params.DemographicParams[agent._race][agent._SO]["PrEPadh"]:
                    agent._PrEP_adh = 1
                else:
                    agent._PrEP_adh = 0
            else:
                if tmp_rnd < params.PrEP_Adherence:
                    agent._PrEP_adh = 1
                else:
                    agent._PrEP_adh = 0

            # set PrEP load and dosestep for PCK
            if params.PrEP_type == "Inj":
                agent._PrEP_load = params.PrEP_peakLoad
                agent._PrEP_lastDose = 0

        if agent == None:
            print("OHHH boi no prep agent")
            return None
        # Check valid input
        if agent._PrEP_bool:
            print((agent._PrEP_bool))
            return None
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

            if params.PrEP_target_model == "Clinical":
                target_PrEP_population = (
                    self.All_agentSet.num_members() - self.HIV_agentSet.num_members()
                )
                target_PrEP = target_PrEP_population * params.PrEP_Target
            elif params.PrEP_target_model == "Incar" or params.PrEP_target_model == "IncarHR":
                if self.runRandom.random() < params.PrEP_Target:
                    _enrollPrEP(self, agent)
                return None
            else:
                target_PrEP = int(
                    (
                        self.All_agentSet.num_members()
                        - self.All_agentSet._subset["HIV"].num_members()
                    )
                    * params.PrEP_Target
                )

            if params.PrEP_clinic_cat == "Racial" and agent_race == "BLACK":
                if self.runRandom.random() < params.PrEP_Target:
                    _enrollPrEP(self, agent)
            elif (
                numPrEP_agents < target_PrEP
                and time >= params.PrEP_startT
                and self._PrEP_elligible(agent, time)
            ):
                # print 'Agent%d added from PrEP'%(agent._ID)
                _enrollPrEP(self, agent)

    def _get_clinic_agent(self, clinicBin, elligiblePool):
        i = 1
        pMatch = params.clinicAgents[clinicBin][i]["Prob"]
        RN = self.runRandom.random()
        while True:
            if RN <= pMatch:
                break
            else:
                i += 1
                pMatch += params.clinicAgents[clinicBin][i]["Prob"]
            if i == 5:
                break

        minNum = params.clinicAgents[clinicBin][i]["min"]
        maxNum = params.clinicAgents[clinicBin][i]["max"]

        iterations = 1
        while iterations < 3:
            randomK_sample = self.runRandom.sample(elligiblePool, params.cal_ptnrSampleDepth)
            elligibleK_Pool = [
                ag
                for ag in randomK_sample
                if ((ag._mean_num_partners >= minNum) and (ag._mean_num_partners <= maxNum))
            ]
            if elligibleK_Pool:
                selected = self.runRandom.choice(elligibleK_Pool)
                elligiblePool.remove(selected)
                return selected
            else:
                print(
                    (
                        "Looking for agent with min:%d and max %d failed %d times"
                        % (minNum, maxNum, iterations)
                    )
                )
                iterations += 1

        print("No suitable PrEP agent")
        return None

    def _progress_to_AIDS(self, agent, agent_drug_type):
        """
        :Purpose:
            Model the progression of HIV agents to AIDS agents
        """
        # only valid for HIV agents
        if not agent._HIV_bool:
            raise ValueError("HAART only valid for HIV agents!agent:%s" % str(agent._ID))

        # if agent not in self.AIDS_agents:
        if not agent._HAART_bool:
            adherenceStat = agent._HAART_adh  # self.AdherenceAgents[agent]
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
                agent._AIDS_bool = True
                self.HIV_AIDS_agentSet.add_agent(agent)


    def _reset_death_count(self):
        self.num_Deaths = {}
        self.deathSet = []
        for HIV_status in ['Total','HIV-', 'HIV+']:
            self.num_Deaths.update({HIV_status: {}})
            for tmp_type in [HIV_status, "MSM", "HM", "HF", "WSW", "MTF"]:
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
        drug_type = self.tmp_Agents[agent]["Drug Type"]

        if drug_type == "IDU":
            try:
                self.tmp_IDU_agents.remove(agent)
            except ValueError:
                pass
        elif drug_type == "NIDU":
            try:
                self.tmp_NIDU_agents.remove(agent)
            except ValueError:
                pass
        elif drug_type == "ND":
            try:
                self.tmp_ND_agents.remove(agent)
            except ValueError:
                pass
        else:
            raise ValueError("Invalid drug type! %s" % str(drug_type))

        # Sex type lists
        try:
            # Agent might have been updated
            sex_type = self.tmp_Agents[agent]["Sex Type"]
        except KeyError:
            sex_type = self.get_agent_characteristic(agent, "Sex Type")
        if sex_type == "MSM":
            self.tmp_MSM_agents.remove(agent)
        elif sex_type == "HF":
            self.tmp_HF_agents.remove(agent)
        elif sex_type == "WSW":
            self.tmp_WSW_agents.remove(agent)
        elif sex_type == "HM":
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
            # print "WTFFFFF"
            pass

        try:
            self.HIVidentified_agents.remove(agent)
        except ValueError:
            # print "WTFFFFF"
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

    def _die_and_replace(self, time, reported=True):
        """
        :Purpose:
            Let agents die and replace the dead agent with a new agent randomly.
        """
        totalDeaths = 0
        for agent in self.All_agentSet._members:

            if agent._incar_bool:
                pass
            else:
                # Probability for dying
                drug_type = agent._DU
                sex_type = agent._SO
                HIV_status = agent._HIV_bool
                HR_status = agent._highrisk_bool
                HR_type = agent._highrisk_type
                AIDSStatus = agent._AIDS_bool
                agent_Race = agent._race

                if HIV_status:
                    if AIDSStatus:  # AIDS DEATH RATE
                        if agent_Race == "WHITE":
                            p = 34.4
                        elif agent_Race == "BLACK":
                            p = 41.6
                        else:
                            raise ValueError("Invalid RACE type! %s" % str(agent_Race))

                    elif agent._HAART_adh > 1:  # HAART DEATH RATE
                        if agent_Race == "WHITE":
                            p = 8.6
                        elif agent_Race == "BLACK":
                            p = 10.4
                        else:
                            raise ValueError("Invalid RACE type! %s" % str(agent_Race))

                    else:  # HIV+ DEATH RATE
                        if agent_Race == "WHITE":
                            p = 17.2
                        elif agent_Race == "BLACK":
                            p = 20.8
                        else:
                            raise ValueError("Invalid RACE type! %s" % str(agent_Race))
                    p = p * params.cal_Mortality

                elif not HIV_status:  # NON HIV DEATH RATE
                    if agent_Race == "WHITE":
                        p = 8.6
                    elif agent_Race == "BLACK":
                        p = 10.4
                    else:
                        raise ValueError("Invalid RACE type! %s" % str(agent_Race))

                else:
                    raise ValueError("Invalid HIV type! %s" % str(HIV_status))

                p = p / 12000.0  #putting it into per 1 person-month #REVIEW does this imply some hardcoded population or timeframe?

                if self.runRandom.random() < p:
                    totalDeaths += 1
                    if HIV_status: ident = "HIV+"
                    else: ident = "HIV-"

                    self.num_Deaths["Total"][sex_type] += 1
                    self.num_Deaths[ident][sex_type] += 1
                    self.deathSet.append(agent)
                    ID_number = agent.get_ID()
                    race = agent._race

                    # End all existing relationships
                    for rel in agent._relationships:
                        rel.progress(forceKill=True)

                        self.Relationships.remove_agent(rel)

                    # Remove agent node and edges from network graph
                    self.get_Graph().remove_node(agent)

                    # Remove agent from agent class and sub-sets
                    self.All_agentSet.remove_agent(agent)

                    # Delete agent object
                    del agent

                    # Create new agent
                    agent_cl = self._return_new_Agent_class(ID_number, race, sex_type)
                    self.create_agent(agent_cl, race)
                    self.G.add_node(agent_cl)

                elif 1 == 0:
                    # Replace with new agent (random characteristics)
                    rv = self.runRandom.random()
                    if rv < params.DemographicParams["WHITE"]["ALL"]["Proportion"]:
                        deliminator = "WHITE"
                    else:
                        deliminator = "BLACK"

                    #################### NOW SET TO REPLAC WITH WHAT DIED"

                    # New agent dict
                    agent_dict = self._return_new_agent_dict(deliminator)
                    if deliminator != agent_dict["Race"]:
                        raise ValueError("Inconsistent drug type!%s" % str(agent_dict["Drug Type"]))

                    # Update tmp_Agents dictionary with new agent
                    self.tmp_Agents.update({agent: agent_dict})

                    # Drug Type
                    drug_type = agent_dict["Drug Type"]
                    if drug_type == "IDU":
                        self.tmp_IDU_agents.append(agent)
                    elif drug_type == "NIDU":
                        self.tmp_NIDU_agents.append(agent)
                    elif drug_type == "ND":
                        self.tmp_ND_agents.append(agent)
                    else:
                        raise ValueError("Invalid drug type! %s" % str(drug_type))

                    # Sex Type
                    SexType = agent_dict["Sex Type"]
                    if SexType == "HM":
                        self.tmp_HM_agents.append(agent)
                    elif SexType == "HF":
                        self.tmp_HF_agents.append(agent)
                    elif SexType == "MSM":
                        self.tmp_MSM_agents.append(agent)
                    elif SexType == "WSW":
                        self.tmp_WSW_agents.append(agent)
                    else:
                        raise ValueError("Invalid SexType! %s" % str(SexType))

                    # HIV
                    HIVStatus = agent_dict["HIV"]
                    if HIVStatus == 1:
                        self.tmp_HIV_agents.append(agent)
                    elif HIVStatus != 0:
                        raise ValueError("Invalid HIVType! %s" % str(HIVStatus))

                    # AIDS
                    AIDSStatus = agent_dict["AIDS"]
                    if AIDSStatus == 1:
                        self.tmp_AIDS_agents.append(agent)
                    elif AIDSStatus != 0:
                        raise ValueError("Invalid AIDS Status! %s" % str(AIDSStatus))

                    # HAART
                    HAARTStatus = agent_dict["HAARTa"]
                    if HAARTStatus == 1:
                        self.tmp_HAART_agents.append(agent)
                    elif HAARTStatus != 0:
                        raise ValueError("Invalid HAART Status! %s" % str(HAARTStatus))

                    # Incarcerated
                    IncarceratedTime = agent_dict["incar_t"]
                    if IncarceratedTime >= 1:
                        self.Incarcerated.append(agent)
                    elif IncarceratedTime < 0:
                        raise ValueError("Invalid AIDS Status! %s" % str(IncarceratedTime))

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

    #REVIEW not used anywhere
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
        self.SEPAgents = {}  # SEP has no memory

    #REVIEW not used anywhere, but may be useful in writing tests?
    def _check_population(self):
        """
        :Purpose:
            Check consistency of population.
            Only called in unittest.

        """

        # Check consistency of last partners
        if not (
            np.all(self.AdjMat.sum(0) == self.AdjMat.conj().sum(0))
            and np.all(self.AdjMat.sum(1) == self.AdjMat.conj().sum(1))
        ):
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
        for (agent, d) in list(self.Agents.items()):
            agent_dict = d
            # Sex type
            sex_type = agent_dict["Sex Type"]
            if sex_type == "HF":
                if agent not in self.HF_agents:
                    print((self.Agents[agent]))
                    raise ValueError("Check agents HF Sex type %d" % agent)
                else:
                    count_HF += 1
            elif sex_type == "HM":
                if agent not in self.HM_agents:
                    print((self.Agents[agent]))
                    raise ValueError("Check agents HM Sex type %d" % agent)
                else:
                    count_HM += 1
            elif sex_type == "MSM":
                if agent not in self.MSM_agents:
                    raise ValueError("Check agents MSM Sex type %d" % agent)
                else:
                    count_MSM += 1
            elif sex_type == "WSW":
                if agent not in self.WSW_agents:
                    print((self.Agents[agent]))
                    raise ValueError("Check agents WSW Sex type %d" % agent)
                else:
                    count_WSW += 1
            else:
                raise ValueError("Invalid sex type %s" % str(sex_type))

            # Drug type
            drug_type = agent_dict["Drug Type"]
            if drug_type == "ND":
                if agent not in self.ND_agents:
                    print((self.Agents[agent]))
                    raise ValueError("Check agents ND Drug type %d" % agent)
                else:
                    count_ND += 1
            elif drug_type == "NIDU":
                if agent not in self.NIDU_agents:
                    print((self.Agents[agent]))
                    raise ValueError("Check agents NIDU Drug type %d" % agent)
                else:
                    count_NIDU += 1
            elif drug_type == "IDU":
                if agent not in self.IDU_agents:
                    print((self.Agents[agent]))
                    raise ValueError("Check agents IDU Drug type %d" % agent)
                else:
                    count_IDU += 1
            else:
                raise ValueError("Invalid drug type %s" % str(drug_type))

            # HIV
            HIVstatus = agent_dict["HIV"]
            if HIVstatus != 0:
                if agent not in self.HIV_agents:
                    print((self.Agents[agent]))
                    raise ValueError("Check agent HIV %d" % agent)
                else:
                    count_HIV += 1
            # AIDS
            AIDSstatus = agent_dict["AIDS"]
            if AIDSstatus != 0:
                if agent not in self.AIDS_agents:
                    print((self.Agents[agent]))
                    raise ValueError("Check agent AIDS %d" % agent)
                else:
                    count_AIDS += 1

        if len(self.HF_agents) != count_HF:
            raise ValueError("self.HF agents contains too many agents!")
        if len(self.HM_agents) != count_HM:
            print(("len(self.HM_agents)=%d" % len(self.HM_agents)))
            print(("count_HM=%d" % count_HM))
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
            raise ValueError(
                "self.HIV_agents contains too many agents!\
                \nlen(self.HIV_agents) = %d\ncount_HIV = %d\n"
                % (len(self.HIV_agents), count_HIV)
            )
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
        for (agent, d) in list(self.tmp_Agents.items()):
            agent_dict = d
            # Sex type
            sex_type = agent_dict["Sex Type"]
            if sex_type == "HF":
                if agent not in self.tmp_HF_agents:
                    print((self.tmp_Agents[agent]))
                    raise ValueError("Check tmp_agents Sex type %d" % agent)
                else:
                    count_HF += 1
            elif sex_type == "HM":
                if agent not in self.tmp_HM_agents:
                    print((self.tmp_Agents[agent]))
                    raise ValueError("Check tmp_agents Sex type %d" % agent)
                else:
                    count_HM += 1
            elif sex_type == "MSM":
                if agent not in self.tmp_MSM_agents:
                    print((self.tmp_Agents[agent]))
                    raise ValueError("Check tmp_agents Sex type %d" % agent)
                else:
                    count_MSM += 1
            elif sex_type == "WSW":
                if agent not in self.tmp_WSW_agents:
                    print((self.tmp_Agents[agent]))
                    raise ValueError("Check tmp_agents Sex type %d" % agent)
                else:
                    count_WSW += 1
            else:
                raise ValueError("Invalid sex type %s" % str(sex_type))

            # Drug type
            drug_type = agent_dict["Drug Type"]
            if drug_type == "ND":
                if agent not in self.tmp_ND_agents:
                    print((self.tmp_Agents[agent]))
                    raise ValueError("Check tmp_agents Drug type %d" % agent)
                else:
                    count_ND += 1
            elif drug_type == "NIDU":
                if agent not in self.tmp_NIDU_agents:
                    print((self.tmp_Agents[agent]))
                    raise ValueError("Check tmp_agents Drug type %d" % agent)
                else:
                    count_NIDU += 1
            elif drug_type == "IDU":
                if agent not in self.tmp_IDU_agents:
                    print((self.tmp_Agents[agent]))
                    raise ValueError("Check tmp_agents Drug type %d" % agent)
                else:
                    count_IDU += 1
            else:
                raise ValueError("Invalid drug type %s" % str(drug_type))

            # HIV
            HIVstatus = agent_dict["HIV"]
            if HIVstatus != 0:
                if agent not in self.tmp_HIV_agents:
                    print((self.tmp_Agents[agent]))
                    raise ValueError("Check tmp_agent HIV %d" % agent)
                else:
                    count_HIV += 1
            # AIDS
            AIDSstatus = agent_dict["AIDS"]
            if AIDSstatus != 0:
                if agent not in self.tmp_AIDS_agents:
                    print((self.tmp_Agents[agent]))
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
            print(("len(self.tmp_AIDS_agents)=%d" % len(self.tmp_AIDS_agents)))
            print(("count_AIDS=%d" % count_AIDS))
            raise ValueError("self.tmp_AIDS agents contains too many agents!")


    #REVIEW not used anywhere
    def save_AgentPartner_list(self, t):
        """
        :Purpsose:
        Save all agent-partners connections.
        :Input:
        t : int
        Time
        """
        OutFileDir = os.path.expanduser(os.path.join(self.current_dir, "Results"))
        if not os.path.isdir(OutFileDir):  # create directory if not existing
            os.mkdir(OutFileDir)
        OutFileName = os.path.join(OutFileDir, "AgentPartnersList_atTime_%s.txt" % str(t))
        if os.path.isfile(OutFileName):
            os.remove(OutFileName)
        outfile = open(OutFileName, "w")
        outfile.write("agent\tdrug type\tsex type\tHIV\tAIDS\tHAART\t")
        maxpartners = 0
        for agent in self.Agents:
            numpartners = len(list(self.AdjMat.rows[agent]))
            if numpartners > maxpartners:
                maxpartners = numpartners
        outfile.write("\t".join(["partner\tp drug type\tp sex type"] * maxpartners))
        outfile.write("\n")
        for agent in sorted(self.Agents.keys()):
            agent_dict = self.Agents[agent]
            outfile.write("%d\t" % agent)
            outfile.write("%s\t" % agent_dict["Drug Type"])
            outfile.write("%s\t" % agent_dict["Sex Type"])
            outfile.write("%d\t" % agent_dict["HIV"])
            outfile.write("%d\t" % agent_dict["AIDS"])
            outfile.write("%d\t" % self.AdherenceAgents[agent])
            for p in sorted(list(self.AdjMat.rows[agent])):
                partner_dict = self.Agents[p]
            outfile.write("%d\t" % int(p))
            outfile.write("%s\t" % partner_dict["Drug Type"])
            outfile.write("%s\t" % partner_dict["Sex Type"])
            outfile.write("\n")


    #REVIEW not used anywhere
    def _reset_partner_count(self):
        """
        Reset partner count for method assess_interaction_distribution
        """
        # set ND partner count to zero for the next time step
        self.tmp_ND_NumPartners_Count = {}
        self.tmp_NIDU_NumPartners_Count = {}
        self.tmp_IDU_NumPartners_Count = {}
        self.tmp_MSM_NumPartners_Count = {}

    #REVIEW not used anywhere
    def get_HIV_prevalence_drugs(self):
        """
        get HIV prevalence within all three drug user groups
        """
        count_HIV_IDU = 0
        count_HIV_NIDU = 0
        count_HIV_ND = 0

        for agent in self.Agents:
            HIVstatus = self.get_agent_characteristic(agent, "HIV")
            if HIVstatus == 1:
                agent_drug_type = self.get_agent_characteristic(agent, "Drug Type")
            if agent_drug_type == "IDU":
                count_HIV_IDU += 1
            elif agent_drug_type == "NIDU":
                count_HIV_NIDU += 1
            elif agent_drug_type == "ND":
                count_HIV_ND += 1
            elif HIVstatus != 0:
                print(HIVstatus)
                raise ValueError("HIV status must be either 0 or 1 !")
                # print [count_HIV_IDU, count_HIV_NIDU, count_HIV_ND]
            else:
                raise ValueError("Agent must be either IDU, NIDU or ND !")
        return [count_HIV_IDU, count_HIV_NIDU, count_HIV_ND]


    #REVIEW not used anywhere
    def get_HIV_prevalence_sex(self):
        """ get HIV prevalence within all four sex groups """
        count_HIV_MSM = 0
        count_HIV_HM = 0
        count_HIV_HF = 0
        count_HIV_WSW = 0

        for agent in self.Agents:
            HIVstatus = self.get_agent_characteristic(agent, "HIV")
            if HIVstatus == 1:
                agent_sex_type = self.get_agent_characteristic(agent, "Sex Type")
            if agent_sex_type == "MSM":
                count_HIV_MSM += 1
            elif agent_sex_type == "HM":
                count_HIV_HM += 1
            elif agent_sex_type == "HF":
                count_HIV_HF += 1
            elif agent_sex_type == "WSW":
                count_HIV_WSW += 1
            elif HIVstatus != 0:
                print(HIVstatus)
                raise ValueError("HIV status must be either 0 or 1 !")
                # print [count_HIV_IDU, count_HIV_NIDU, count_HIV_ND]
            else:
                raise ValueError("Agent must be either MSM, HM, MF, or WSW !")

        return [count_HIV_MSM, count_HIV_HM, count_HIV_HF, count_HIV_WSW]

    #REVIEW not used anywhere
    def get_HIV_prevalence_drugs_sex(self):
        """prevalences without and msm only"""
        count_HIV_MIDU = 0
        count_HIV_MNIDU = 0
        count_HIV_MND = 0
        count_HIV_IDUnmsm = 0
        count_HIV_NIDUnmsm = 0
        count_HIV_NDnmsm = 0

        for agent in self.Agents:
            HIVstatus = self.get_agent_characteristic(agent, "HIV")
            if HIVstatus == 1:
                agent_sex_type = self.get_agent_characteristic(agent, "Sex Type")
            agent_drug_type = self.get_agent_characteristic(agent, "Drug Type")
            if agent_drug_type == "IDU" and agent_sex_type in ["HM", "HF", "WSW"]:
                count_HIV_IDUnmsm += 1
            elif agent_drug_type == "IDU" and agent_sex_type == "MSM":
                count_HIV_MIDU += 1
            elif agent_drug_type == "NIDU" and agent_sex_type in ["HM", "HF", "WSW"]:
                count_HIV_NIDUnmsm += 1
            elif agent_drug_type == "NIDU" and agent_sex_type == "MSM":
                count_HIV_MNIDU += 1
            elif agent_drug_type == "ND" and agent_sex_type in ["HM", "HF", "WSW"]:
                count_HIV_NDnmsm += 1
            elif agent_drug_type == "ND" and agent_sex_type == "MSM":
                count_HIV_MND += 1
            elif HIVstatus != 0:
                print(HIVstatus)
            raise ValueError("HIV status must be either 0 or 1 !")
        return [
            count_HIV_MIDU,
            count_HIV_MNIDU,
            count_HIV_MND,
            count_HIV_IDUnmsm,
            count_HIV_NIDUnmsm,
            count_HIV_NDnmsm,
        ]


    #REVIEW not used anywhere
    def get_HIV_prevalence(self):
        """ get HIV prevalence"""
        HIVcount = 0.0
        for agent in list(self.Agents.keys()):
            HIVstatus = self.get_agent_characteristic(agent, "HIV")
            if HIVstatus == 1:
                HIVcount += 1
        return HIVcount


    def return_results(self):
        return self.ResultDict

    #REVIEW not used anywhere
    def save_result_dict(self):
        OutFileDir = os.path.join(self.current_dir, "Results")
        if not os.path.isdir(OutFileDir):  # create directory if not existing
            os.mkdir(OutFileDir)
        OutFileName = os.path.join(OutFileDir, "ResultDictionary.txt")
        if os.path.isfile(OutFileName):
            os.remove(OutFileName)
        outfile = open(OutFileName, "w")
        for result_property in sorted(self.ResultDict.keys()):
            outfile.write("%s\t" % result_property)
            for time_t in sorted(self.ResultDict[result_property].keys()):
                outfile.write("%4.5f\t" % float(self.ResultDict[result_property][time_t]))
            outfile.write("\n")


    #REVIEW not used anywhere
    def save_AdjMat(self, t):
        """
        :Purpose:
        Save Adjacency matrix in sparse format.
        """
        OutFileDir = os.path.expanduser(os.path.join(self.current_dir, "Results"))
        if not os.path.isdir(OutFileDir):  # create directory if not existing
            os.mkdir(OutFileDir)
        OutFileName = os.path.join(OutFileDir, "AdjacencyMatrix_atTime_%d.txt" % t)
        if os.path.isfile(OutFileName):
            os.remove(OutFileName)
        outfile = open(OutFileName, "w")
        for n, row in enumerate(self.AdjMat.rows):
            outfile.write("%d:\t" % n)
            for partner in row:
                outfile.write("%d," % partner)
            outfile.write("\n")

    #REVIEW not used anywhere
    def _writeDNR(self):
        dynnetworkReport = open("Results/dynnetworkReport.txt", "a")
        for agent in self.Agents:
            sextype = self.Agents[agent]["Sex Type"]
            drugtype = self.Agents[agent]["Drug Type"]
            HIV = self.Agents[agent]["HIV"]
            reportLine = "\t".join(["0", "NEWAGENT", repr(agent), sextype, drugtype, repr(HIV)])
            dynnetworkReport.write("\n" + reportLine)

        dynnetworkReport.close()
        # initiate HIV status
        for agent in self.Agents:
            if self.Agents[agent]["HIV"] == 1:
                self._become_HIV(agent, 0)
