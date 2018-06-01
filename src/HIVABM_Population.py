#!/usr/bin/env python2.3
"""
*****************************************************************************
Author(s):	Maximilian King  (previous authors: Lars Seemann - lseemann@uh.edu)
Email: Maximilian_King@brown.edu
Organization: Marshall Lab, Department of Epidemiology - Brown University


Description:
    Module for population dictionary construction and initial parameters.
    Contains population demographics for model inputs, and creates respective
    population subsets.

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

import random
import copy
from copy import deepcopy
import unittest

try:
    from agent import *
except ImportError:
    raise ImportError("Can't import AgentClass")
import params

def print_population(agent_dict, dir_prefix, time=0):
    # Print the whole population to a file.
    # Format: agent_i characteristic 1 characteristic 2 ...
    # Input: agent_dict : dictionary of agents (self.Agents)
    #        dir_prefix : outfile name (str)
    #        time : time at which snapshot was taken (int)
    OutFileName = (dir_prefix + '/' +'Population_at_time_%d.txt'%time)
    outfile = open(OutFileName,'w')
    # create header
    key_list = agent_dict[0].keys()
    outfile.write('Agent\t')
    for key in key_list:
        outfile.write('%s\t'%str(key))
    outfile.write('\n')
    for agent in sorted(agent_dict.keys()):
        outfile.write('%d\t'%agent)
        tmp_agent_prop = agent_dict[agent]
        for key in key_list:
            outfile.write('%s\t'%str(tmp_agent_prop[key]))
        outfile.write('\n')
    outfile.close()

class PopulationClass():

    """
    :Purpose:
        This class constructs and represents the model population

    :Input:	
    
        n : int
            Number of agents. Default: 10000

    :Attributes:
    
        :py:attr:`PopulationSize` : int
            Size of the population.

        :py:attr:`propIDU` : float
            Prevalence of intravenous drug users.

        :py:attr:`numIDU` : int
            Number of intravenous drug users.

        :py:attr:`propNIDU` : float
            Prevalence of non-intravenous drug users.

        :py:attr:`numNIDU` : int
            Number of non-intravenous drug users.

        :py:attr:`propND` : float
            Prevalence of not drug users.

        :py:attr:`numND` : int
            Number of not drug users.

        :py:attr:`Agents`: dict
            Dictionary of agents and their characteristics. The agents are
            the `keys' and a dicotionary of 'characteristic:value'
            pair is the entry.

        :py:attr:`IDU_agents`: list
            IDU drug users

        :py:attr:`NIDU_agents`: list
            NIDU drug users

        :py:attr:`ND_agents`: list
            ND drug users

        :py:attr:`MSM_agents`: list
            MSM agents

        :py:attr:`HIV_agents`: list
            HIV+ users

        :py:attr:`AIDS_agents`: list
            Users with AIDS

    :Methods:
        :py:meth:`get_agent_characteristic` \n
        :py:meth:`_set_population` \n
        :py:meth:`get_agents` \n
        :py:meth:`get_info_DrugUserType` \n
        :py:meth:`get_info_HIV_IDU` \n
        :py:meth:`get_info_DrugSexType`
    """

    def __init__(self, n = 10000, rSeed = 0, model=None):
        """
        :Purpose:
            Initialize PopulationClass object.
        """
        self.RANDOMSEED = rSeed
        print rSeed
        #random.seed(self.RANDOMSEED)
        if type(n) is not int:
            raise ValueError("Population size must be integer")
        else:
            self.PopulationSize = n
        # Parameters

        #Stratification proportions
        self.propWhite = 1
        self.propBlack = 0


        self.numWhite = int(params.DemographicParams['WHITE']['ALL']['Proportion'] * self.PopulationSize)
        self.numBlack = int(params.DemographicParams['BLACK']['ALL']['Proportion'] * self.PopulationSize)
        print 'W: %d, B: %d'%(self.numWhite, self.numBlack)
        #Nested dictionary for probability lookups by race
        #StratW = White, StratB = Black
        #HM incar @ 0.0279


        #MODEL TYPES
        MT_PrEP = False
        MT_Incar = False
        MT_NoIncar = False

        if model == 'PrEP':
            MT_PrEP = True
        elif model == 'Incar':
            MT_Incar = True
        elif model == 'NoIncar':
            MT_NoIncar = True

        #Quick calc based on traj: 0.596 Female -> 0,01175/0.75 HIV -> 0.75 Tested+ = 11.75% HIV DIAG
        StratW = {'MSM':{}, 'HM':{}, 'HF':{}, 'PWID':{}}
        if MT_PrEP:
            StratW['MSM'] = {'POP':1.0, 'HIV':0.036195675, 'AIDS':0.51, 'HAARTprev':0.33, 'INCARprev':0.0, 'TestedPrev':0.816, 'mNPart':3}
            StratW['HM'] = {'POP':0.0, 'HIV':0.0, 'AIDS':0.0, 'HAARTprev':0.0, 'INCARprev':0.0, 'TestedPrev':0.0}
            StratW['HF'] = {'POP':0.0, 'HIV':0.0, 'AIDS':0.0, 'HAARTprev':0.0, 'INCARprev':0.0, 'TestedPrev':0.0}
            StratW['PWID'] = {'POP':0.0, 'HIV':0.0, 'AIDS':0.0, 'HAARTprev':0.0, 'INCARprev':0.0, 'TestedPrev':0.0}
        elif MT_Incar:
            StratW['MSM'] = {'POP':0.0, 'HIV':0.05, 'AIDS':0.05, 'HAARTprev':0.0, 'INCARprev':0.0, 'TestedPrev':0.86, 'mNPart':3}
            StratW['HM'] = {'POP':0.4044435412, 'HIV':0.0158/0.75, 'AIDS':0.678, 'HAARTprev':0.03748, 'INCARprev':0.0279, 'TestedPrev':0.75, 'mNPart':3}
            StratW['HF'] = {'POP':0.5955564588, 'HIV':0.01175/0.75, 'AIDS':0.573, 'HAARTprev':0.04054, 'INCARprev':0.0, 'TestedPrev':0.75, 'mNPart':2}
            StratW['PWID'] = {'POP':0.0, 'HIV':0.0, 'AIDS':0.0, 'HAARTprev':0.0, 'INCARprev':0.0, 'TestedPrev':0.0, 'mNPart':3}
        elif MT_NoIncar:
            StratW['MSM'] = {'POP':0.0, 'HIV':0.05, 'AIDS':0.05, 'HAARTprev':0.0, 'INCARprev':0.0, 'TestedPrev':0.86, 'mNPart':3}
            StratW['HM'] = {'POP':0.4044435412, 'HIV':0.0158/0.75, 'AIDS':0.678, 'HAARTprev':0.03748, 'INCARprev':0.0, 'TestedPrev':0.75, 'mNPart':3}
            StratW['HF'] = {'POP':0.5955564588, 'HIV':0.01175/0.75, 'AIDS':0.573, 'HAARTprev':0.04054, 'INCARprev':0.0, 'TestedPrev':0.75, 'mNPart':2}
            StratW['PWID'] = {'POP':0.0, 'HIV':0.0, 'AIDS':0.0, 'HAARTprev':0.0, 'INCARprev':0.0, 'TestedPrev':0.0, 'mNPart':3}

        StratB = {'MSM':{}, 'HM':{}, 'HF':{}, 'PWID':{}}
        StratB['MSM'] = {'POP':0.0, 'HIV':0.0 , 'AIDS':0.0, 'HAARTprev':0.0, 'INCARprev':0.0, 'TestedPrev':0.0}
        StratB['HM'] = {'POP':0.0, 'HIV':0.0 , 'AIDS':0.0, 'HAARTprev':0.0, 'INCARprev':0.0, 'TestedPrev':0.0}
        StratB['HF'] = {'POP':0.0, 'HIV':0.0 , 'AIDS':0.0, 'HAARTprev':0.0, 'INCARprev':0.0, 'TestedPrev':0.0}
        StratB['PWID'] = {'POP':0.0, 'HIV':0.0 , 'AIDS':0.0, 'HAARTprev':0.0, 'INCARprev':0.0, 'TestedPrev':0.0}

        self.ProbLookUp = {'WHITE':StratW, 'BLACK':StratB}
        # drug user prevalence (proportion)


        self.propIDU  = 0##190/10000.0
        self.numIDU = int(self.propIDU * self.PopulationSize)
        self.propNIDU = 0##640/10000.0
        self.numNIDU = int(self.propNIDU * self.PopulationSize)
        self.propND = 0.995##9170/10000.0
        self.numND = self.PopulationSize - self.numIDU - self.numNIDU
        self.Agents = dict()	# Main Dictionary Agents
        # List of IDU, NIDU, NDs
        # shuffle all Agents
        allAgents = range(self.PopulationSize)
        #random.shuffle(allAgents)

        print "\tBuilding class sets"
        self.totalAgentClass = Agent_set(0,"TotalAgents")
        self.MSM_agentsClass = Agent_set(1,"MSM", numerator=self.totalAgentClass)
        self.HM_agentsClass = Agent_set(1,"HM", numerator=self.totalAgentClass)
        self.HF_agentsClass = Agent_set(1,"HF", numerator=self.totalAgentClass)
        self.MTF_agentsClass = Agent_set(1,"MTF", numerator=self.totalAgentClass)
        self.IDU_agentsClass = Agent_set(1,"IDU", numerator=self.totalAgentClass)
        self.HIV_agents_class = Agent_set(2,"HIV", numerator=self.totalAgentClass)
        self.HAART_agentClass = Agent_set(2,"HAART", numerator=self.HIV_agents_class)
        self.PrEP_agents_class = Agent_set(2,"PrEP")
        self.IncarceratedClass = Agent_set(2,"Incar", numerator=self.totalAgentClass)
        self.HighriskClass = Agent_set(2,"HRisk", numerator=self.totalAgentClass)

        print "\tOrganizing subsets and heirarchy"
        self.totalAgentClass.add_subset(self.HIV_agents_class)
        self.totalAgentClass.add_subset(self.PrEP_agents_class)
        self.totalAgentClass.add_subset(self.IncarceratedClass)
        self.totalAgentClass.add_subset(self.HAART_agentClass)
        self.totalAgentClass.add_subset(self.HighriskClass)
        self.totalAgentClass.add_subset(self.MSM_agentsClass)
        self.totalAgentClass.add_subset(self.HM_agentsClass)
        self.totalAgentClass.add_subset(self.HF_agentsClass)
        self.totalAgentClass.add_subset(self.MTF_agentsClass)
        self.totalAgentClass.add_subset(self.IDU_agentsClass)


        self.Relationships = Agent_set(1,"Relationships")
        """### TESTING NEW CLAS
        agent1 = Agent(1, "HM", 28, "White", "ND")

        self.totalAgents = Agent_set(2,1)
        self.totalAgents.add_agent(agent1)
        self.totalAgents.add_agent(Agent(2,"HF",26,"White","ND"))
        self.totalAgents.add_agent(Agent(3,"HF",26,"White","ND"))
        print self.totalAgents.is_member(2)
        print self.totalAgents.is_member(5)
        print self.totalAgents.iter_agents()

        self.testAgents = Agent_set(2,1)
        for tmp in self.totalAgents.iter_agents():
            if tmp._SO == "HM":self.testAgents.add_agent(tmp)

        self.testAgents.print_agents()
        self.totalAgents.print_agents()

        temp = self.totalAgents.get_agent(1)
        temp._SO = "MSM"

        self.testAgents.print_agents()
        self.totalAgents.print_agents()
        """



        # Nested dictionary for probability look-up
        IDU = {'MSM':{}, 'HM':{}, 'WSW':{}, 'HF':{}}
        IDU['MSM'] = {'HIV':0.55 , 'AIDS':0.058,'HAARTa':0}
        IDU['HM'] = {'HIV':0.42 , 'AIDS':0.058, 'HAARTa':0}
        IDU['WSW'] = {'HIV':0.53 , 'AIDS':0.058, 'HAARTa':0}
        IDU['HF'] = {'HIV':0.39 , 'AIDS':0.058, 'HAARTa':0}
        NIDU = {'MSM':{}, 'HM':{}, 'WSW':{}, 'HF':{}}
        NIDU['MSM'] = {'HIV':0.18 , 'AIDS':0.02, 'HAARTa':0}
        NIDU['HM'] = {'HIV':0.048 , 'AIDS':0.002, 'HAARTa':0}
        NIDU['WSW'] = {'HIV':0.048 , 'AIDS':0.002, 'HAARTa':0}
        NIDU['HF'] = {'HIV':0.048 , 'AIDS':0.002, 'HAARTa':0}
        ND = {'MSM':{}, 'HM':{}, 'WSW':{}, 'HF':{}}
        ND['Type'] = ([0,1,2,3], [0.469, 0.493, 0.022, 0.016])
        ND['MSM'] = {'HIV':0.08 , 'AIDS':0.02, 'HAARTa':0}
        ND['HM'] = {'HIV':0.015 , 'AIDS':0.0003, 'HAARTa':0}
        ND['WSW'] = {'HIV':0.012 , 'AIDS':0.0003, 'HAARTa':0}
        ND['HF'] = {'HIV':0.012 , 'AIDS':0.0003, 'HAARTa':0}
        #self.ProbLookUp = {'IDU':IDU, 'NIDU':NIDU1, 'ND':ND}

        # print ('Number of agents '+str(len(allAgents)))


        self.MSM_agents = []
        self.HF_agents = []
        self.WSW_agents = []
        self.HM_agents = []
        self.IDU_agents = []
        self.HIV_agents = []
        self.AIDS_agents = []
        self.HIVidentified_agents = []
        self.HAART_agents = []
        self.Incarcerated = []

        #Create agents in allAgents list
        #"""
        self.White_agents = deepcopy(allAgents[0:self.numWhite])
        self.Black_agents = deepcopy(allAgents[self.numWhite:])
        self.IDU_agents = []#deepcopy(allAgents[0:self.numIDU])
        self.NIDU_agents = []#deepcopy(allAgents[self.numIDU:(self.numIDU + self.numNIDU)])
        self.ND_agents = []#deepcopy(allAgents[(self.numIDU + self.numNIDU):])

        print "\tCreating agents"

        for agent in self.White_agents:
            self.create_agent(agent, 'WHITE')
        for agent in self.Black_agents:
            self.create_agent(agent, 'BLACK')
        random.shuffle(self.totalAgentClass._members)
        prob_Incarc = params.DemographicParams['WHITE']['HM']['INCARprev']
        for tmpA in self.totalAgentClass._members:

            dice = random.random()
            prob_Incarc = params.DemographicParams[tmpA._race][tmpA._SO]['INCARprev']
            if dice < prob_Incarc:
                #print dice, prob_Incarc
                toss = 2#random.choice( (1, 2) )
                if toss == 1: #JAIL
                    tmpA._incar_bool = True
                    tmpA._incar_time = int(random.triangular(1, 9, 3))
                else: #PRISON
                    tmpA._incar_bool = True
                    tmpA._incar_time = int(random.triangular(1, params.inc_PrisMax, int(params.inc_PrisMax/3))) #random.randint(params.inc_PrisMin, params.inc_PrisMax)
                self.IncarceratedClass.add_agent(tmpA)
            #self.totalAgentClass._subset["Incar"].add_agent(agent_cl)
        #"""
        #self.totalAgentClass.print_agents()
        # [self.MSM_agentsClass.add_agent(tmpA) for tmpA in self.totalAgentClass.iter_agents() if tmpA._SO == "MSM"]
        # [self.HM_agentsClass.add_agent(tmpA) for tmpA in self.totalAgentClass.iter_agents() if tmpA._SO == "HM"]
        # [self.HF_agentsClass.add_agent(tmpA) for tmpA in self.totalAgentClass.iter_agents() if tmpA._SO == "HF"]
        """for tmpA in self.totalAgentClass.iter_agents():
            if tmpA._SO == "HM":
                self.MSM_agentsClass.add_agent(tmpA)
            elif tmpA._SO == "MSM":
                self.MSM_agentsClass.add_agent(tmpA)
            elif tmpA._SO == "HF":
                self.MSM_agentsClass.add_agent(tmpA)
            else:
                print "EERRROROOOR"
                exit()
            if tmpA._HIV_bool:
                self.totalAgentClass._subset["HIV"].add_agent(tmpA)
            """



        """
        for agent in self.Agents:
            racetype = self.Agents[agent]['Race']
            sextype = self.Agents[agent]['Sex Type']
            drugtype = self.Agents[agent]['Drug Type']
            print ('%.6d\t28\tU\t%s\t%s\t%s'%(agent, sextype,drugtype, racetype))



        self.IDU_agents = deepcopy(allAgents[0:self.numIDU])
        self.NIDU_agents = deepcopy(allAgents[self.numIDU:(self.numIDU + self.numNIDU)])
        self.ND_agents = deepcopy(allAgents[(self.numIDU + self.numNIDU):])

        for agent in self.IDU_agents:
            self.create_agent(agent,DrugType = 'IDU')
        for agent in self.NIDU_agents:
            self.create_agent(agent,DrugType = 'NIDU')
        for agent in self.ND_agents:
            self.create_agent(agent,DrugType = 'ND')

        # Check consistency
        CheckSum_SexType = len(self.MSM_agents)+len(self.HF_agents)+len(self.WSW_agents)+len(self.HM_agents)
        if CheckSum_SexType != self.PopulationSize:
            raise ValueError("MSM:%d\nHF:%d\nWSW:%d\nHM:%d\nSum:%d"%(
                len(self.MSM_agents),len(self.HF_agents),len(self.WSW_agents),len(self.HM_agents),CheckSum_SexType))

        CheckSum_DrugType = len(self.IDU_agents)+len(self.NIDU_agents)+len(self.ND_agents)
        if CheckSum_DrugType != self.PopulationSize:
            raise ValueError("IDU:%d\nNIDU:%d\nND:%d\nSum:%d"%(
                len(self.IDU_agents),len(self.NIDU_agents),len(self.ND_agents),CheckSum_DrugType))
        """

    def _return_agent_set(self):
        return self.totalAgentClass


    def get_agent_characteristic(self, agent, property):
        """
        :Purpose:
            Determine an agents interal characteristics.
        :Input:
            agent : int \n
            property: string
        	Allowed properties: \n
        	`HIV`, `Sex Type`, `Drug Type`
        :Output:
            Characteristic : string
        """
        if property in ['HIV', 'Sex Type', 'Drug Type', 'AIDS', 'HAARTa', 'incar_t', 'HIV_t', 'Race', 'Tested']:
            try:
                return self.Agents[agent][property]
            except KeyError:
                print self.Agents[agent]
                raise KeyError("Check keys/Agents in self.Agents")
        else:
            raise ValueError("Check agents charactertic! Only 'HIV', 'Sex Type' \
			     'Drug Type' , 'AIDS' possible.")

    def _return_new_agent_dict(self, Deliminator):
        """
        :Purpose:
        Return random agent dict of a new agent..
            Each agent is a key to an associated dictionary which stores the internal
            characteristics in form of an additinoal dictionary of the form
            ``characteristic:value``.
        :Input:
        DrugType : str
            Either 'IDU','NIDU' or 'ND'
        :Output:
             agent_dict : dict
        """
        SexType = 'NULL'
        Drugtype = 'NULL'

        #Determine sextype
        tmp_rnd = random.random()
        if tmp_rnd < params.DemographicParams[Deliminator]['HM']['POP']:
            SexType = 'HM'
        elif tmp_rnd < (params.DemographicParams[Deliminator]['HM']['POP'] + params.DemographicParams[Deliminator]['HF']['POP']):
            SexType = 'HF'
        else:
            SexType = 'MSM'

        #Determine drugtype
        tmp_rnd = random.random()
        #print "%.3lf must be less than%.3lf"%(tmp_rnd,params.DemographicParams[Deliminator]['PWID']['POP'])

        if tmp_rnd < params.DemographicParams[Deliminator]['PWID']['POP']:
            DrugType = 'IDU'
        else:
            DrugType = 'ND'

        # HIV
        if DrugType == 'IDU':
            prob_HIV = params.DemographicParams[Deliminator]['PWID']['HIV']
        else:
            prob_HIV = params.DemographicParams[Deliminator][SexType]['HIV']
        #print "%s\t%s\t%.3lf"%(Deliminator,SexType,prob_HIV)

        if random.random() < prob_HIV:
            HIVStatus = 1

            # if HIV AIDS possible
            if DrugType == 'IDU':
                prob_AIDS = params.DemographicParams[Deliminator]['PWID']['AIDS']
            else:
                prob_AIDS = params.DemographicParams[Deliminator][SexType]['AIDS']
                #print "%s\t%s\t%.3lf"%(Deliminator,SexType,prob_HIV)

            if random.random() < prob_AIDS:
                AIDSStatus = 1
            else:
                AIDSStatus = 0

            # HIV testing params
            if DrugType == 'IDU':
                prob_Tested = params.DemographicParams[Deliminator]['PWID']['TestedPrev']
            else:
                prob_Tested = params.DemographicParams[Deliminator][SexType]['TestedPrev']

            if random.random() < prob_Tested:
                TestedStatus = 1

                #if tested HAART possible
                if DrugType == 'IDU':
                    prob_HAART = params.DemographicParams[Deliminator]['PWID']['HAARTprev']
                else:
                    prob_HAART = params.DemographicParams[Deliminator][SexType]['HAARTprev']

                if random.random() < prob_HAART:
                    HAARTStatus = 1
                else:
                    HAARTStatus = 0

            else:
                TestedStatus = 0
                HAARTStatus = 0

            # if HIV, how long has the agent had it? Random sample
            HIV_time = random.randint(1,42)

        else:
            HIVStatus = 0
            AIDSStatus = 0
            HAARTStatus = 0
            HIV_time = 0
            TestedStatus = 0

        #Incarceration
        if DrugType == 'IDU':
            prob_Incarc = params.DemographicParams[Deliminator]['PWID']['INCARprev']
        else:
            prob_Incarc = params.DemographicParams[Deliminator][SexType]['INCARprev']

        if random.random() < prob_Incarc:
            toss = random.choice( (1, 2) )
            if toss == 1: #JAIL
                incar_time = int(random.triangular(6, 21, 15))
            else: #PRISON
                incar_time = int(random.triangular(30, 96, 60))
        else:
            incar_time = 0

        #print "New agent: %s\t%s\t%s\tHIV:%d" % (Deliminator,DrugType,SexType,HIVStatus)
        agent_dict = {'Race':Deliminator,'Drug Type': DrugType,'Sex Type':SexType, 'HIV':HIVStatus, 'Tested':TestedStatus, 'AIDS':AIDSStatus, 'HAARTa':HAARTStatus, 'incar_t':incar_time,'HIV_t':HIV_time}

        return agent_dict

    def _return_new_Agent_class(self, agentID, Race):
        """
        :Purpose:
        Return random agent dict of a new agent..
            Each agent is a key to an associated dictionary which stores the internal
            characteristics in form of an additinoal dictionary of the form
            ``characteristic:value``.
        :Input:
        DrugType : str
            Either 'IDU','NIDU' or 'ND'
        :Output:
             agent_dict : dict
        """
        SexType = 'NULL'
        Drugtype = 'NULL'

        #Determine sextype
        tmp_rnd = random.random()
        if tmp_rnd < params.DemographicParams[Race]['HM']['POP']:
            SexType = 'HM'
        elif tmp_rnd < (params.DemographicParams[Race]['HM']['POP'] + params.DemographicParams[Race]['HF']['POP']):
            SexType = 'HF'
        elif tmp_rnd < (params.DemographicParams[Race]['HM']['POP'] + params.DemographicParams[Race]['HF']['POP'] + params.DemographicParams[Race]['MSM']['POP']):
            SexType = 'MSM'
        else:
            SexType = 'MTF'

        #Determine drugtype
        tmp_rnd = random.random()
        #print "%.3lf must be less than%.3lf"%(tmp_rnd,params.DemographicParams[Deliminator]['PWID']['POP'])

        #todo: FIX THIS TO GET BACK IDU
        if tmp_rnd < params.DemographicParams[Race]['PWID']['POP']:
            DrugType = 'IDU'
        else:
            DrugType = 'ND'


        #age = random.randint(18,65)

        age, ageBin = self.getAge(Race)


        newAgent = Agent(agentID, SexType, age, Race, DrugType)
        newAgent._ageBin = ageBin
        # HIV
        if DrugType == 'IDU':
            prob_HIV = params.DemographicParams[Race]['PWID']['HIV']
        else:
            #prob_HIV = params.ageMatrix[Race]['HIV'][newAgent._ageBin]
            prob_HIV = params.DemographicParams[Race][SexType]['HIV']
        #print "%s\t%s\t%.3lf"%(Drugtype,SexType,prob_HIV)

        if random.random() < prob_HIV:
            HIVStatus = 1
            newAgent._HIV_bool = True

            # if HIV AIDS possible
            if DrugType == 'IDU':
                prob_AIDS = params.DemographicParams[Race]['PWID']['AIDS']
            else:
                prob_AIDS = params.DemographicParams[Race][SexType]['AIDS']
                #print "%s\t%s\t%.3lf"%(Deliminator,SexType,prob_HIV)

            if random.random() < prob_AIDS:
                AIDSStatus = 1
                newAgent._AIDS_bool = True
            else:
                AIDSStatus = 0

            # HIV testing params
            if DrugType == 'IDU':
                prob_Tested = params.DemographicParams[Race]['PWID']['TestedPrev']
            else:
                prob_Tested = params.DemographicParams[Race][SexType]['TestedPrev']

            if random.random() < prob_Tested:
                TestedStatus = 1
                newAgent._tested = True

                #if tested HAART possible
                if DrugType == 'IDU':
                    prob_HAART = params.DemographicParams[Race]['PWID']['HAARTprev']
                else:
                    prob_HAART = params.DemographicParams[Race][SexType]['HAARTprev']

                if random.random() < prob_HAART:
                    HAARTStatus = 1
                    newAgent._HAART_bool = True
                else:
                    HAARTStatus = 0

            else:
                TestedStatus = 0
                HAARTStatus = 0

            # if HIV, how long has the agent had it? Random sample
            #HIV_time = random.randint(1,42)
            newAgent._HIV_time = random.randint(1,42)

        else:

            HIVStatus = 0
            AIDSStatus = 0
            HAARTStatus = 0
            HIV_time = 0
            TestedStatus = 0

        # #Incarceration
        # if DrugType == 'IDU':
        #     prob_Incarc = params.DemographicParams[Race]['PWID']['INCARprev']
        # else:
        #     prob_Incarc = params.DemographicParams[Race][SexType]['INCARprev']
        #
        # if random.random() < prob_Incarc:
        #     toss = random.choice( (1, 2) )
        #     if toss == 1: #JAIL
        #         newAgent._incar_bool = True
        #         newAgent._incar_time = int(random.triangular(1, 9, 3))
        #     else: #PRISON
        #         newAgent._incar_bool = True
        #         newAgent._incar_time = int(random.triangular(6, 60, 24))
        # else:
        #     incar_time = 0


        diceroll = random.random()

        if diceroll < 0.01:
            mNPart = 0
        elif diceroll < (0.01 + 0.14):
            mNPart = 1
        elif diceroll < (0.01 + 0.14 + 0.16):
            mNPart = 2
        elif diceroll < (0.01 + 0.14 + 0.16 + 0.19):
            mNPart = random.randrange(3,4,1)
        elif diceroll < (0.01 + 0.14 + 0.16 + 0.19 + 0.34):
            mNPart = random.randrange(5,10,1)
        else: #16%
            mNPart = 10
        #Partnership demographics
        newAgent._mean_num_partners = mNPart #params.DemographicParams[Race][SexType]['mNPart']

        #print "New agent: %s\t%s\t%s\tHIV:%d" % (Deliminator,DrugType,SexType,HIVStatus)
        #agent_dict = {'Race':Race,'Drug Type': DrugType,'Sex Type':SexType, 'HIV':HIVStatus, 'Tested':TestedStatus, 'AIDS':AIDSStatus, 'HAARTa':HAARTStatus, 'incar_t':incar_time,'HIV_t':HIV_time}

        return newAgent

    def create_agent(self, agent, Deliminator):
        """
	:Purpose:
	    Creat a new agent in the population.
            Each agent is a key to an associated dictionary which stores the internal 
            characteristics in form of an additinoal dictionary of the form
            ``characteristic:value``.
	
	:Input:
	    agent : int
	    
	    DrugType : str
	        Either 'IDU','NIDU' or 'ND'
	"""

        agent_cl = self._return_new_Agent_class(agent,Deliminator)
        self.totalAgentClass.add_agent(agent_cl)
        if agent == 0.25*self.PopulationSize:print "25%"
        elif agent == 0.5*self.PopulationSize:print "50%"
        elif agent == 0.75*self.PopulationSize:print "75%"
        #elif agent == self.PopulationSize:print "100%"

        if agent_cl._DU == 'IDU':
            self.totalAgentClass._subset["IDU"].add_agent(agent_cl)

        if agent_cl._HIV_bool:
            self.totalAgentClass._subset["HIV"].add_agent(agent_cl)

        if agent_cl._incar_bool:
            self.IncarceratedClass.add_agent(agent_cl)
            #self.totalAgentClass._subset["Incar"].add_agent(agent_cl)

        if agent_cl._HAART_bool:
            self.HAART_agentClass.add_agent(agent_cl)


        """
        agent_dict = self._return_new_agent_dict(Deliminator)

        if Deliminator != agent_dict['Race']:
            raise ValueError("Inconsistent drug type!%s"%str(
                agent_dict['Race']))

        # Agent dict
        self.Agents.update({agent:agent_dict})

        # Drug Type
        DrugType = agent_dict['Drug Type']
        if DrugType == 'IDU':
            self.IDU_agents.append(agent)
        elif DrugType == 'ND':
            self.ND_agents.append(agent)
        else:
            raise ValueError("Invalid SexType! %s"%str(DrugType))

        # Sex Type
        SexType = agent_dict['Sex Type']
        if SexType == 'HM':
            self.HM_agents.append(agent)
        elif SexType == 'HF':
            self.HF_agents.append(agent)
        elif SexType == 'MSM':
            self.MSM_agents.append(agent)
        elif SexType == 'WSW':
            self.WSW_agents.append(agent)
        else:
            raise ValueError("Invalid SexType! %s"%str(SexType))

        # HIV
        HIVStatus = agent_dict['HIV']
        if HIVStatus == 1:
            self.HIV_agents.append(agent)
        elif HIVStatus != 0:
            raise ValueError("Invalid HIVType! %s"%str(HIVStatus))

        # AIDS
        AIDSStatus = agent_dict['AIDS']
        if AIDSStatus == 1:
            self.AIDS_agents.append(agent)
        elif AIDSStatus != 0:
            raise ValueError("Invalid AIDS Status! %s"%str(AIDSStatus))

        # HAART
        HAARTStatus = agent_dict['HAARTa']
        if HAARTStatus == 1:
            self.HAART_agents.append(agent)
        elif HAARTStatus != 0:
            raise ValueError("Invalid HAART Status! %s"%str(HAARTStatus))

        #Incarcerated
        IncarceratedTime = agent_dict['incar_t']
        if IncarceratedTime >= 1:
            self.Incarcerated.append(agent)
        elif IncarceratedTime < 0:
            raise ValueError("Invalid AIDS Status! %s"%str(IncarceratedTime))

        # TESTED
        TestStatus = agent_dict['Tested']
        if TestStatus == 1:
            self.HIVidentified_agents.append(agent)
        elif TestStatus != 0:
            raise ValueError("Invalid HAART Status! %s"%str(HAARTStatus))
        """

    def getAge(self, race):
        rand = random.random()
        minAge = 15
        maxAge = 80
        ageBin = 0
        #print params.ageMatrix[race]['Prop'][1]

        if rand < params.ageMatrix[race]['Prop'][1]:
            minAge = 15
            maxAge = 24
            ageBin = 1
        elif rand < params.ageMatrix[race]['Prop'][2]:
            minAge = 25
            maxAge = 34
            ageBin = 2
        elif rand < params.ageMatrix[race]['Prop'][3]:
            minAge = 35
            maxAge = 44
            ageBin = 3
        elif rand < params.ageMatrix[race]['Prop'][4]:
            minAge = 45
            maxAge = 54
            ageBin = 4
        else:
            minAge = 55
            maxAge = 80
            ageBin = 5
        # else:
        #     minAge = 15
        #     maxAge = 80
        age = random.randrange(minAge,maxAge)
        #print rand, race, age, ageBin
        return age, ageBin

    def get_agents(self):
        """ 
        :Purpose:
            Return all agents and their characteristics.

        :Output: 
            :py:attr:`Agents`: dict
        	Dictionary of agents and their characteristics. The agents
        	are the `keys` and a dictionary of `characteristic:value`
        	pair is the entry. 
        """
        return self.Agents

    def print_info(self):
        """ 
        :Purpose:
            Simple fprintf test on std out.
        """
        print ('Number of agents '+str(self.PopulationSize))
        print ('Number of IDU '+str(len(self.IDU_agents)))
        print ('Number of NIDU '+str(len(self.NIDU_agents)))
        print ('Number of ND '+str(len(self.ND_agents)))
        """
	# random IDU
	agent = random.choice(self.IDU.keys())
	print (' Random IDU agent: ' + str(agent))
	print self.IDU[agent]
	# random NIDU
	agent = random.choice(self.NIDU.keys())
	print (' Random NIDU agent: ' + str(agent))
	print self.NIDU[agent]
	# random ND
	agent = random.choice(self.ND.keys())
	print (' Random ND agent: ' + str(agent))
	print self.ND[agent]
	"""

    def get_info_DrugUserType(self):
        """ 
        :Purpose:
            Return number of IDU, NIDU, and ND users in one array.

        :Output:
            data : dict
        """
        return {'Number of IDU agents': len(self.IDU_agents),
                'Number of NIDU agents':len(self.NIDU_agents),
                'Number of ND agents':len(self.ND_agents)}

    def get_info_HIV_IDU(self):
        """
        :Purpose: 
            Return number of HIV among IDU agents. 
            Distinguish between MSM, WSW, HM, and HIF agents.

        :Ooutput: 
            data : dict
        """
        count_HIV_MSM = 0
        count_HIV_HM  = 0
        count_HIV_WSW = 0
        count_HIV_HF  = 0
        for agent in self.IDU_agents:
            HIVstatus = self.get_agent_characteristic(agent, 'HIV')
            if HIVstatus == 1:
                SexType = self.get_agent_characteristic(agent, 'Sex Type')
                if SexType == 'MSM': count_HIV_MSM += 1
                elif SexType == 'HM':count_HIV_HM += 1
                elif SexType == 'WSW':count_HIV_WSW += 1
                elif SexType == 'HF':count_HIV_HF += 1
                else:raise ValueError("Agent must be either HM, MSM, HF, WSW !")
            elif HIVstatus != 0:
                print HIVstatus
                raise ValueError("HIV status must be either 0 or 1 !")
        return {'Number of IDU HIV MSM':count_HIV_MSM,
                'Number of IDU HIV HM':count_HIV_HM,
                'Number of IDU HIV WSW':count_HIV_WSW,
                'Number of IDU HIV HF':count_HIV_HF}

    def get_info_DrugSexType(self):
        """
        :Purpose: 
            Assess the Drug and Sex type prevalences of the population.

        :Ooutput: 
            data : dict
        """
        count_HM	= 0
        count_HF	= 0
        count_MSM	= 0
        count_WSW	= 0
        count_HM_IDU	= 0
        count_HF_IDU	= 0
        count_MSM_IDU	= 0
        count_WSW_IDU	= 0
        count_HM_NIDU	= 0
        count_HF_NIDU	= 0
        count_MSM_NIDU	= 0
        count_WSW_NIDU	= 0
        count_HM_ND	= 0
        count_HF_ND	= 0
        count_MSM_ND	= 0
        count_WSW_ND	= 0

        for agent in self.Agents.keys():
            SexType = self.get_agent_characteristic(agent, 'Sex Type')
            DrugType = self.get_agent_characteristic(agent, 'Drug Type')
            if SexType == 'HM':
                count_HM += 1
                if DrugType == 'IDU':count_HM_IDU += 1
                elif DrugType == 'NIDU':count_HM_NIDU += 1
                elif DrugType == 'ND':count_HM_ND += 1
                else: raise ValueError("Drug test must either be IDU, NIDU or ND !")
            elif SexType == 'HF':
                count_HF += 1
                if DrugType == 'IDU': count_HF_IDU += 1
                elif DrugType == 'NIDU':count_HF_NIDU += 1
                elif DrugType == 'ND':count_HF_ND += 1
                else:raise ValueError("Drug test must either be IDU, NIDU or ND !")
            elif SexType == 'MSM':
                count_MSM += 1
                if DrugType == 'IDU':count_MSM_IDU += 1
                elif DrugType == 'NIDU':count_MSM_NIDU += 1
                elif DrugType == 'ND':count_MSM_ND += 1
                else:raise ValueError("Drug test must either be IDU, NIDU or ND !")
            elif SexType == 'WSW':
                count_WSW += 1
                if DrugType == 'IDU':count_WSW_IDU += 1
                elif DrugType == 'NIDU':count_WSW_NIDU += 1
                elif DrugType == 'ND':count_WSW_ND += 1
                else:raise ValueError("Drug test must either be IDU, NIDU or ND !")
            else:raise ValueError("Agent must be either HM, MSM, HF, WSW !")

        return {'Number of HM':count_HM,
                'Number of HF':count_HF,
                'Number of MSM':count_MSM,
                'Number of WSW':count_WSW,
                'Number of MSM IDU':count_MSM_IDU,
                'Number of MSM NIDU':count_MSM_NIDU,
                'Number of MSM ND':count_MSM_ND,
                'Number of HM IDU':count_HM_IDU,
                'Number of HM NIDU':count_HM_NIDU,
                'Number of HM ND':count_HM_ND,
                'Number of WSW IDU':count_WSW_IDU,
                'Number of WSW NIDU':count_WSW_NIDU,
                'Number of WSW ND':count_WSW_ND,
                'Number of HF IDU':count_HF_IDU,
                'Number of HF NIDU':count_HF_NIDU,
                'Number of HF ND':count_HF_ND}

class TestClassMethods(unittest.TestCase):
    """ 
    :Purpose:
        unittest
    """
    def setUp(self):
        """
        :Purpose:
            Tests that all models from setup pass inspection. ``setUp`` is perfomed before each method.
        """
        self.N_pop = 10000

    def test_SexType(self):
        """ Test: Testing consistency of Sex type agents"""
        print "\t__Testing HM agents list"
        myPopulation = PopulationClass(n=self.N_pop)
        for a in myPopulation.Agents.keys():
            SexType = myPopulation.get_agent_characteristic(a,'Sex Type')
            if SexType == 'HM':
                self.assertTrue(a in myPopulation.HM_agents)
            elif SexType == 'HF':
                self.assertTrue(a in myPopulation.HF_agents)
            elif SexType == 'MSM':
                self.assertTrue(a in myPopulation.MSM_agents)
            elif SexType == 'WSW':
                self.assertTrue(a in myPopulation.WSW_agents)

        for agent in myPopulation.HM_agents:
            self.assertTrue(myPopulation.get_agent_characteristic(agent,'Sex Type') == 'HM')
        for agent in myPopulation.HF_agents:
            self.assertTrue(myPopulation.get_agent_characteristic(agent,'Sex Type') == 'HF')
        for agent in myPopulation.MSM_agents:
            self.assertTrue(myPopulation.get_agent_characteristic(agent,'Sex Type') == 'MSM')
        for agent in myPopulation.WSW_agents:
            self.assertTrue(myPopulation.get_agent_characteristic(agent,'Sex Type') == 'WSW')

    def test_HIV(self):
        """ Test: Testing HIV agent array"""
        print "\t__Testing the HIV agent list"
        tmpCount = 0
        myPopulation = PopulationClass(n=self.N_pop)
        for a in myPopulation.Agents.keys():
            HIVstatus = myPopulation.get_agent_characteristic(a,'HIV')
            if HIVstatus != 0:
                tmpCount += 1
                self.assertTrue(a in myPopulation.HIV_agents)
        self.assertTrue(len(myPopulation.HIV_agents) == tmpCount)

    def test_AIDS(self):
        """ Test: Testing AIDS agent array"""
        print "\t__Testing the AIDS agent list"
        tmpCount = 0
        myPopulation = PopulationClass(n=self.N_pop)
        for a in myPopulation.Agents.keys():
            AIDSstatus = myPopulation.get_agent_characteristic(a,'AIDS')
            if AIDSstatus != 0:
                tmpCount += 1
                self.assertTrue(a in myPopulation.AIDS_agents)
                self.assertTrue(a in myPopulation.HIV_agents)
        self.assertTrue(len(myPopulation.AIDS_agents) == tmpCount)

    def test_consistency(self):
        """ Test: Testing consistency"""
        print "\t__Testing consistency of agent lists"
        myPop = PopulationClass(n=self.N_pop)
        NormalAgents = list(set(range(myPop.PopulationSize)).difference(
            set(myPop.IDU_agents).union(set(myPop.MSM_agents))))
        MSMAgents = list(set(myPop.MSM_agents).difference(set(myPop.IDU_agents)))
        IDUagents = myPop.IDU_agents
        NumAgents = len(NormalAgents)+len(MSMAgents)+len(IDUagents)
        self.assertTrue(NumAgents == self.N_pop,"NumAgents=%d, \
				PopulationSize = %d"%(NumAgents, self.N_pop))

    def test_MSM(self):
        """ Test: Testing MSM agent array"""
        print "\t__Testing the MSM agent list"
        tmpCount = 0
        myPopulation = PopulationClass(n=self.N_pop)
        for a in myPopulation.Agents.keys():
            SEXstatus = myPopulation.get_agent_characteristic(a,'Sex Type')
            if SEXstatus == 'MSM':
                tmpCount += 1
                self.assertTrue(a in myPopulation.MSM_agents)
        self.assertTrue(len(myPopulation.MSM_agents) == tmpCount)

        for agent in myPopulation.MSM_agents:
            agent_sex_type = myPopulation.get_agent_characteristic(agent,'Sex Type')
            self.assertTrue(agent_sex_type == 'MSM')

    def test_IDU(self):
        """ Test: Testing IDU agent array"""
        print "\t__Testing the IDU agent list"
        tmpCount = 0
        myPopulation = PopulationClass(n=self.N_pop)
        for agent in myPopulation.Agents.keys():
            agent_drug_status = myPopulation.get_agent_characteristic(agent,'Drug Type')
            if agent_drug_status == 'IDU':
                tmpCount += 1
                self.assertTrue(agent in myPopulation.IDU_agents)
        self.assertTrue(len(myPopulation.IDU_agents) == tmpCount)

        for agent in myPopulation.IDU_agents:
            agent_drug_type = myPopulation.get_agent_characteristic(agent,'Drug Type')
            self.assertTrue(agent_drug_type == 'IDU')

    def test_NIDU(self):
        """ Test: Testing INDU agent array"""
        print "\t__Testing the NIDU agent list"
        tmpCount = 0
        myPopulation = PopulationClass(n=self.N_pop)
        for agent in myPopulation.Agents.keys():
            agent_drug_status = myPopulation.get_agent_characteristic(agent,'Drug Type')
            if agent_drug_status == 'NIDU':
                tmpCount += 1
                self.assertTrue(agent in myPopulation.NIDU_agents)
        self.assertTrue(len(myPopulation.NIDU_agents) == tmpCount)

        for agent in myPopulation.NIDU_agents:
            agent_drug_type = myPopulation.get_agent_characteristic(agent,'Drug Type')
            self.assertTrue(agent_drug_type == 'NIDU')

    def test_Population(self):
        """ Test: Testing the population"""
        print "\t__Testing the population"
        myPopulation = PopulationClass(n=self.N_pop)
        for agent in myPopulation.Agents.keys():
            tmp_DrugType = myPopulation.get_agent_characteristic(agent,'Drug Type')
            self.assertTrue(tmp_DrugType in ['NIDU','IDU','ND'])
            if tmp_DrugType == 'NIDU':
                self.assertTrue(agent in myPopulation.NIDU_agents)
            elif tmp_DrugType == 'IDU':
                self.assertTrue(agent in myPopulation.IDU_agents)
            else: self.assertTrue(agent in myPopulation.ND_agents)
        self.assertTrue(len(myPopulation.NIDU_agents) == myPopulation.numNIDU)
        self.assertTrue(len(myPopulation.IDU_agents) == myPopulation.numIDU)
        self.assertTrue(len(myPopulation.ND_agents) == myPopulation.numND)

if __name__=='__main__':
    """unittest"""
    unittest.main()				
