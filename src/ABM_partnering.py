#!/usr/bin/env python
# encoding: utf-8

"""
*****************************************************************************
Author(s):	Maximilian King  (previous authors: Lars Seemann - lseemann@uh.edu)
Email: Maximilian_King@brown.edu
Organization: Marshall Lab, Department of Epidemiology - Brown University

Description:
    Module responsible for partnering agents within network. Assortative mixing,
    partnering preferences, and eligible lists.


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
from copy import deepcopy, copy
import os
import time
#import PyQt4

#import scipy.sparse as spsp
from scipy.stats import binom
from scipy.stats import poisson
from functools import wraps
import numpy as np

try:
    from HIVABM_Population import PopulationClass, print_population
except ImportError:
    raise ImportError("Can't import PopulationClass")

try:
    from ABM_core import *
except ImportError:
    raise ImportError("Can't import PopulationClass")

try:
    from agent import *
except ImportError:
    raise ImportError("Can't import Agent class")

import params

def update_partner_assignments(self, partnerTurnover, graph, agent=None):
    # Now create partnerships until available partnerships are out

    if agent:
        partner = get_partner(self, agent, self.totalAgentClass)

        if partner:
            #print "Agent %d found partner %d!"%(agent.get_ID(), partner.get_ID())
            duration = get_partnership_duration(self, agent)
            tmp_relationship = Relationship(agent, partner, "MSM", "SE", duration)
            agent.bond(partner, tmp_relationship)
            self.Relationships.add_agent(tmp_relationship)
            graph.G.add_edge(tmp_relationship._ID1, tmp_relationship._ID2)
    else:
        EligibleAgents = self.totalAgentClass#._subset["HIV"].iter_agents()
        noMatch = 0

        for agent in EligibleAgents.iter_agents():
            #print len(agent._partners)
            acquirePartnerProb = (agent._mean_num_partners / (12.0))*partnerTurnover
            #if agent._highrisk_bool:print acquirePartnerProb
            if np.random.uniform(0, 1) < acquirePartnerProb:
                partner = get_partner(self, agent, self.totalAgentClass)

                if partner:
                    #print "Agent %d found partner %d!"%(agent.get_ID(), partner.get_ID())

                    duration = get_partnership_duration(self, agent)
                    tmp_relationship = Relationship(agent, partner, "MSM", "SE", duration)
                    agent.bond(partner, tmp_relationship)
                    self.Relationships.add_agent(tmp_relationship)
                    graph.G.add_edge(tmp_relationship._ID1, tmp_relationship._ID2)
                    # ADD RELATIONSHIP EDGE TO GRAPH G of NetworkGraph
                    #print "%d/%d partnets found for agent %d"%(len(agent._partners), agent._num_sex_partners, agent.get_ID())
                    #print "%d/%d partnets found for partner %d"%(len(partner._partners), partner._num_sex_partners, partner.get_ID())

                else:
                    #print "Missed pass attempt",noMatch
                    noMatch += 1

    # print "\n\t\t-COULDNT MATCH",noMatch,"AGENTS IN NEED \t---"

def get_number_of_partners(self, agent, agent_drug_type, agent_sex_type):
    """
    :Purpose:
        Get number of partners for a agent.
        Drawn from Poisson distribution.

    :Input:
        agent_drug_type : str
        Either 'IDU', 'NIDU', 'ND'

        agent_sex_type : str
        Either 'HM', 'MSM', 'HF', 'WSW'

    :Output:
        NumPartners : int
        Zero partners possible.
    """
    # Check input
    # Drug type
    if agent_drug_type not in ['IDU', 'NIDU', 'ND']:
        raise ValueError("Invalid drug type! %s" % str(agent_drug_type))
    # Sex type
    if agent_sex_type not in ['HM', 'HF', 'MSM', 'WSW']:
        raise ValueError("Invalid sex type! %s" % str(agent_sex_type))

    agent_race_type = agent._race

    n_trials = self.ProbTables[agent_race_type][agent_sex_type]['NUMPartn']
    p_success = .8

    ##Random number of contacts using negative binomial
    if agent_sex_type == 'WSW':
        # n_trials = 1
        # p_success = 0.8
        RandNumCont = np.random.negative_binomial(n_trials, p_success, 1)[0]
    elif agent_sex_type == 'MSM' and agent_drug_type != 'NIDU':
        # n_trials = 1
        # p_success = 0.8
        RandNumCont = np.random.negative_binomial(n_trials, p_success, 1)[0]
    elif agent_sex_type == 'MSM' and agent_drug_type == 'NIDU':
        # n_trials = 1
        # p_success = 0.8
        RandNumCont = np.random.negative_binomial(n_trials, p_success, 1)[0]
        RandNumCont = int(RandNumCont * 2)
    elif agent_drug_type == 'NIDU':
        # n_trials = 1
        # p_success = 0.8
        RandNumCont = np.random.negative_binomial(n_trials, p_success, 1)[0]
    elif agent_drug_type == 'IDU':
        n_trials = 7
        p_success = 0.7
        RandNumCont = np.random.negative_binomial(n_trials, p_success, 1)[0]
    elif agent_drug_type == 'ND':
        # n_trials = 1
        # p_success= 0.8
        RandNumCont = np.random.negative_binomial(n_trials, p_success, 1)[0]
    if RandNumCont < 0:
        raise ValueError("Invalid number of contacts!%s" % str(RandNumCont))

    if RandNumCont == 0 and np.random.uniform() < .5:
        RandNumCont = 1

    MEAN_PARTNER_YEAR = self.ProbTables[agent_race_type][agent_sex_type]['NUMPartn']
    RandNumCont = poisson.rvs(MEAN_PARTNER_YEAR, size=1)

    if agent_drug_type == 'IDU':
        RandNumCont = RandNumCont * params.cal_NeedlePartScaling

    else:
        RandNumCont = RandNumCont * params.cal_SexualActScaling


    return RandNumCont


def get_partner(self, agent, need_new_partners):
    """
    :Purpose:
        Get partner for agent.

    :Input:
        agent : int

        need_new_partners: list of available partners

    :Output:
        partner: new partner
    """
    #print need_new_partners
    shortlist_NNP = need_new_partners

    agent_sex_type = agent._SO
    agent_drug_type = agent._DU
    RandomPartner = None


    #print("Finding partner for agent", agent._ID, agent_sex_type, agent_drug_type)
    if agent_drug_type == 'IDU':
        if random.random() < 0.8:
            # choose from IDU agents
            try: RandomPartner = get_random_IDU_partner(self, agent, shortlist_NNP)
            except:
                print "No IDU matches"
                get_random_sex_partner(self, agent, shortlist_NNP)
            #print "\tReturned: %s" % RandomPartner
        else:
            get_random_sex_partner(self, agent, shortlist_NNP)
    elif agent_drug_type in ('ND','NIDU'):
        if params.flag_AssortativeMix:
            if random.random() < params.AssortMixCoeff:
                RandomPartner = get_assort_sex_partner(self, agent, shortlist_NNP)
            else:
                RandomPartner = get_random_sex_partner(self, agent, shortlist_NNP)
        else:
            RandomPartner = get_random_sex_partner(self, agent, shortlist_NNP)
    else:
        raise ValueError("Check method _get_partners(). Agent not caught!")
    #print RandomPartner

    #RandomPartner = random.choice(need_new_partners)
    if RandomPartner == agent: return None
    else: return RandomPartner


def get_random_IDU_partner(self, agent, need_new_partners):
    """
    :Purpose:
        Get a random partner which is sex compatible

    :Input:
        agent: int
        need_new_partners: list of available partners

    :Output:
        partner : int

    """
    agent_sex_type = agent._SO
    agent_drug_type = agent._DU
    RandomPartner = None
    tempList = []

    AssortMix = False
    if random.random() < params.AssortMixCoeff:
        AssortMix = True
    #assert agent_drug_type in ['IDU'], "Invalid drug type for IDU! %s"%str(agent_drug_type)
    #todo: Make the random agent never return the agent or any of their partners
    if agent_drug_type not in ['IDU']:
        raise ValueError("Invalid drug type! %s"%str(agent_drug_type))
    else:    
        RandomPartner = random.choice(need_new_partners._subset["IDU"]._members)
        if RandomPartner in agent._partners or RandomPartner == agent:
            RandomPartner = None

    #print "\tReturned: %s" % RandomPartner
    if RandomPartner:
        return RandomPartner
    else:
        return None
        #print "NO PATNEAS"

def get_assort_IDU_partner(self, agent, need_new_partners, assortType):
    """
    :Purpose:
        Get a random partner which is sex compatible and adherese to assortativity mixing defined by params

    :Input:
        agent: int
        need_new_partners: list of available partners

    :Output:
        partner : int

    """
    agent_sex_type = agent._SO
    agent_drug_type = agent._DU
    RandomPartner = None
    tempList = []

    AssortMix = False
    if random.random() < params.AssortMixCoeff:
        AssortMix = True
    #assert agent_drug_type in ['IDU'], "Invalid drug type for IDU! %s"%str(agent_drug_type)
    #todo: Make the random agent never return the agent or any of their partners
    if agent_drug_type not in ['IDU']:
        raise ValueError("Invalid drug type! %s"%str(agent_drug_type))
    else:
        RandomPartner = random.choice(need_new_partners._subset["IDU"]._members)
        if RandomPartner in agent._partners or RandomPartner == agent:
            RandomPartner = None

    #print "\tReturned: %s" % RandomPartner
    if RandomPartner:
        return RandomPartner
    else:
        return None
        #print "NO PATNEAS"

def get_assort_sex_partner(self, agent, need_new_partners):
    """
    :Purpose:
        Get a random partner which is sex compatible and fits assortativity constraints

    :Input:
        agent: int
        need_new_partners: list of available partners

    :Output:
        partner : int

    """
    def partner_choice(x):
        intersection = list(set(need_new_partners).intersection(set(x)))
        agent_race_type = self.get_agent_characteristic(agent, 'Race')
        #print agent_race_type
        if agent_race_type == 'WHITE':
            Assortive_intersection = list(set(self.White_agents).intersection(intersection))
            if Assortive_intersection == []: print "Couldnt assortive mix (W), picking suitable agent"
            else: return random.choice(Assortive_intersection)
        elif agent_race_type == 'BLACK':
            Assortive_intersection = list(set(self.Black_agents).intersection(intersection))
            if Assortive_intersection == []:
                print "Couldnt assortive mix (B), picking suitable agent"
            else:
                #print Assortive_intersection
                return random.choice(Assortive_intersection)
        if intersection == []: return None
        else: print "NO PATNAS"#return random.choice(intersection)

    def getPartnerBin(agent):

        testRand = random.random()
        i = 1
        pMatch = params.mixingMatrix[agent._ageBin][i]

        #print params.mixingMatrix[1][1]
        #print agent._ageBin, i
        if params.flag_AgeAssortMix:
            while(True):
                if testRand <= pMatch:
                    return i
                else:
                    i+=1
                    pMatch += params.mixingMatrix[agent._ageBin][i]
                if i==5:return i
        else:
            i = random.randrange(1,6)
            return i


    agent_sex_type = agent._SO
    agent_drug_type = agent._DU
    agent_race_type = agent._race

    RandomPartner = None
    tempList = []

    if random.random() < params.AssortMixCoeff:
        AssortMix = True
    else:
        AssortMix = False

    rv = random.random()

    #todo: Make the random agent never return the agent or any of their partners
    assert(agent_sex_type in ['HM','HF','MSM','WSW','MTF'])

    eligPartnerType = params.DemographicParams[agent_race_type][agent_sex_type]['EligPartnerType'][0]

    if params.AssortMixType == 'Age':
        randomK_sample = random.sample(need_new_partners._subset["MSM"]._members, params.cal_ptnrSampleDepth)
        ageBinPick = getPartnerBin(agent)
        while True:
            RandomPartner = random.choice([ag for ag in randomK_sample if ag._ageBin == ageBinPick])
            break

    #else if picking using race mix
    elif params.AssortMixType == 'Race':
        samplePop = [tmpA for tmpA in need_new_partners._subset[eligPartnerType]._members if tmpA._race != agent._race]
        try:
            randomK_sample = random.sample(samplePop ,params.cal_ptnrSampleDepth)
        except:
            randomK_sample = samplePop
        while True:
            RandomPartner = random.choice([ag for ag in randomK_sample if ag._race == agent._race])
            break

    elif params.AssortMixType == 'HR':
        samplePop = [tmpA for tmpA in need_new_partners._subset[eligPartnerType]._members if tmpA._everhighrisk_bool]
        if samplePop:
            try:
                randomK_sample = random.sample(samplePop ,params.cal_ptnrSampleDepth)
            except:
                randomK_sample = samplePop

            while True:
                RandomPartner = random.choice([ag for ag in randomK_sample if ag._race == agent._race])
                break

    if RandomPartner == None or RandomPartner in agent._partners or RandomPartner == agent:
        RandomPartner = None
    else:
        try: RandomPartner = random.choice(tempList)
        except: pass#print "No matches in", tempList

    #print "\tReturned: %s" % RandomPartner
    if RandomPartner:
        return RandomPartner
    else:
        pass
        #print "NO PATNEAS"

def get_random_sex_partner(self, agent, need_new_partners):
    """
    :Purpose:
        Get a random partner which is sex compatible

    :Input:
        agent: int
        need_new_partners: list of available partners

    :Output:
        partner : int

    """
    def partner_choice(x):
        intersection = list(set(need_new_partners).intersection(set(x)))
        agent_race_type = self.get_agent_characteristic(agent, 'Race')
        #print agent_race_type
        if agent_race_type == 'WHITE':
            Assortive_intersection = list(set(self.White_agents).intersection(intersection))
            if Assortive_intersection == []: print "Couldnt assortive mix (W), picking suitable agent"
            else: return random.choice(Assortive_intersection)
        elif agent_race_type == 'BLACK':
            Assortive_intersection = list(set(self.Black_agents).intersection(intersection))
            if Assortive_intersection == []:
                print "Couldnt assortive mix (B), picking suitable agent"
            else:
                #print Assortive_intersection
                return random.choice(Assortive_intersection)
        if intersection == []: return None
        else: print "NO PATNAS"#return random.choice(intersection)

    def getPartnerBin(agent):

        testRand = random.random()
        i = 1
        pMatch = params.mixingMatrix[agent._ageBin][i]

        #print params.mixingMatrix[1][1]
        #print agent._ageBin, i
        if params.flag_AgeAssortMix:
            while(True):
                if testRand <= pMatch:
                    return i
                else:
                    i+=1
                    pMatch += params.mixingMatrix[agent._ageBin][i]
                if i==5:return i
        else:
            i = random.randrange(1,6)
            return i


    #agent_sex_type = self.get_agent_characteristic(agent, 'Sex Type')
    agent_sex_type = agent._SO
    agent_race_type = agent._race
    agent_drug_type = agent._DU
    #print "\tChecking for sex partner for %d" % agent
    RandomPartner = None
    tempList = []

    eligPartnerType = params.DemographicParams[agent_race_type][agent_sex_type]['EligPartnerType'][0]
    AssortMix = False
    if params.flag_AgeAssortMix:
        if random.random() < params.AssortMixCoeff:
            AssortMix = True

    #todo: Make the random agent never return the agent or any of their partners
    if agent_sex_type not in ['HM','HF','MSM','WSW','MTF']:
        raise ValueError("Invalid sex type! %s"%str(agent_sex_type))
    else:
        while True:
            RandomPartner = random.choice(need_new_partners._subset[eligPartnerType]._members)

            if RandomPartner in agent._partners or RandomPartner == agent:
                pass
            else:
                break
        if RandomPartner in agent._partners or RandomPartner == agent:
            RandomPartner = None


    #print "\tReturned: %s" % RandomPartner
    if RandomPartner:
        assert(sex_possible(self,agent._SO, RandomPartner._SO)),"Sex no possible between agents! ERROR 441"
        return RandomPartner
    else:
        pass


def sex_possible(self, agent_sex_type, partner_sex_type):
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
    if agent_sex_type not in ['HM', 'HF', 'MSM', 'WSW','MTF']:
        raise ValueError("Invalid agent_sex_type! %s" % str(agent_sex_type))
    if partner_sex_type not in ['HM', 'HF', 'MSM', 'WSW','MTF']:
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


def get_partnership_duration(self, agent):
    """
    :Purpose:
        Get number of partners for a agent.
        Drawn from Poisson distribution.

    :Input:
        agent_drug_type : agentClass

    :Output:
        NumPartners : int
        Zero partners possible.
    """
    # Check input
    agent_drug_type = agent._DU
    agent_sex_type = agent._SO

    # Drug type
    if agent_drug_type not in ['IDU', 'NIDU', 'ND']:
        raise ValueError("Invalid drug type! %s" % str(agent_drug_type))
    # Sex type
    if agent_sex_type not in params.agentSexTypes:
        raise ValueError("Invalid sex type! %s" % str(agent_sex_type))

    diceroll = random.random()

     # Length of relationship (months)a
    # <1 1,679 32.3% 566 17.7 1,113 55.8
    # 1–6 1,359 26.2% 929 29.0 430 21.6
    # 7–12 604 11.6% 459 14.4 145 7.3
    # 13–24 628 12.1% 480 15.0 148 7.4
    # 25–36 309 6.0% 264 8.3 45 2.3
    # >37 614 11.8% 501 15.7 113 5.7

    if diceroll < params.sexualDurations[1]['p_value']:
        dur_bin = 1
    elif diceroll < params.sexualDurations[2]['p_value']:
        dur_bin = 2
    elif diceroll < params.sexualDurations[3]['p_value']:
        dur_bin = 3
    elif diceroll < params.sexualDurations[4]['p_value']:
        dur_bin = 4
    else:
        dur_bin = 5

    duration = random.randrange(params.sexualDurations[dur_bin]['min'], params.sexualDurations[dur_bin]['max'], 1)

    return duration

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



def update_partner_assignment_old(self, partnerTurnover):
    # Generate target partner numbers for each agent and get current partner nums
    target_partner_nums = {}
    current_partner_nums = {}

    print "ZZZZZZZZZZZZ"

    # Loop through agents again, if too few: go into need_partners set
    need_new_partners = []
    need_new_partners2 = Agent_set(2,"NNP")

    need_new_partners_HM = Agent_set(2,"HM")
    need_new_partners_MSM = Agent_set(2,"MSM")
    need_new_partners_HF = Agent_set(2,"HF")

    need_new_partners2.add_subset(need_new_partners_HM)
    need_new_partners2.add_subset(need_new_partners_HF)
    need_new_partners2.add_subset(need_new_partners_MSM)

    """
    # TODO: FIX THIS BACK TO HIV AGENTS ONLY
    EligibleAgents = self.totalAgentClass.iter_agents()
    for agent in EligibleAgents:
        agent_sex_type = agent._SO
        agent_drug_type = agent._DU

        current_num = len(agent._partners)
        #print current_num

        if np.random.uniform(0, 1) > partnerTurnover:
            target_num = current_num
        else:
            target_num = get_number_of_partners(self, agent, agent_drug_type, agent_sex_type)
            agent._num_sex_partners = target_num



        #p_success = 12/(12+agent._mean_num_partners)
        #newPartnerProb = np.random.negative_binomial(12, p_success)
        #print newPartnerProb

        # Now loop through agents, if currently too many partners, remove some
        if target_num < current_num:
            # ExistingLinks = list(self.AdjMat.rows[agent])
            n = current_num - target_num
            for i in range(n):
                agent2remove = random.choice(agent._partners)
                # print "Current agent %d has %d partners and wants %d - " %(agent, current_partner_nums[agent],target_partner_nums[agent]), list(self.AdjMat.rows[agent]), "removing %d"%agent2remove
                # self.AdjMat[agent, agent2remove] = 0  # remove connection adjMat
                # self.AdjMat[agent2remove, agent] = 0
                # current_partner_nums[agent] -= 1
                # if agent2remove in EligibleAgents:  # self.HIV_agents:
                #     current_partner_nums[agent2remove] -= 1
                #print "Removing %d from %s (agent %d)"%(agent2remove.get_ID(), agent.partner_list(), agent.get_ID())
                agent.unpair(agent2remove)




    # need_new_partners2._subset["HM"] = need_new_partners_HM
    # need_new_partners2._subset["MSM"] = need_new_partners_MSM
    # need_new_partners2._subset["HF"] = need_new_partners_HF
    #print need_new_partners2._subset["HM"]

    #todo: fix this to only use the neednewpartner list
    EligibleAgents = self.totalAgentClass.iter_agents()
    for agent in EligibleAgents:#self.HIV_agents:#Agents:
        #print len(agent._partners)

        if agent._num_sex_partners > len(agent._partners):
            # need_new_partners.append(agent)
            need_new_partners2.add_agent(agent)
            #print agent._SO
            need_new_partners2._subset[agent._SO].add_agent(agent)


            # if agent._SO == "HM":
            #     need_new_partners2._subset["HM"].add_agent(agent)
            # elif agent._SO == "MSM":
            #     need_new_partners2._subset["MSM"].add_agent(agent)
            # elif agent._SO == "HF":
            #     need_new_partners2._subset["HF"].add_agent(agent)

            #print agent
            #print need_new_partners

    #need_new_partners = list(np.random.permutation(need_new_partners))
    #need_new_partners2.print_agents()
    """

    # Now create partnerships until available partnerships are out
    EligibleAgents = self.totalAgentClass.iter_agents()
    self.Relationships
    for agent in EligibleAgents:#self.HIV_agents:#Agents:
        #print len(agent._partners)
        acquirePartnerProb = agent._mean_num_partners / 12.0
        #print acquirePartnerProb
        if np.random.uniform(0, 1) < acquirePartnerProb:
            need_new_partners2.add_agent(agent)
            #print agent._SO
            need_new_partners2._subset[agent._SO].add_agent(agent)
    print "\t\t-FINDING MATCHES FOR",need_new_partners2.num_members(),"AGENTS IN NEED \t---"

    #need_new_partners2._subset["HM"].print_agents()
    noMatch = 0
    missed_counter=0
    print_counter = need_new_partners2.num_members()
    print_percent = 1.0

    while need_new_partners2.num_members() > 0:

        #need_new_partners2.print_agents()
        #need_new_partners2.print_subsets()
        # if need_new_partners2.num_members() <= int(print_counter * print_percent):
        #     print "Remain %d" % (print_percent*100),
        #     print_percent -= 0.1
        # else:
        #     print "Dafuq"
        #print "NNP Remaining:",need_new_partners2.num_members()
        agent = None
        #self.totalAgentClass.print_agents()
        #print need_new_partners2.print_agents()
        #print "Selected agent ..."
        agent = need_new_partners2.random_agent()
        #agent.print_agent()
        # TODO Fix this to read proper agent partner

        partner = None
        while partner == None:
            #print "\nGetting partner for agent %d"%agent.get_ID()
            partner = get_partner(self, agent, need_new_partners2)
            #partner_cl = self.totalAgentClass.get_agent(partner)
            #self.AdjMat[agent, partner] = 1
            #self.AdjMat[partner, agent] = 1
            #current_partner_nums[agent] += 1

            if partner:
                #print "Agent %d found partner %d!"%(agent.get_ID(), partner.get_ID())

                duration = 10
                agent.bond(partner, duration)
                tmp_relationship = Relationship(agent, partner, "MSM", "SE", duration)
                self.Relationships.add_agent(tmp_relationship)
                #print "%d/%d partnets found for agent %d"%(len(agent._partners), agent._num_sex_partners, agent.get_ID())
                #print "%d/%d partnets found for partner %d"%(len(partner._partners), partner._num_sex_partners, partner.get_ID())
                #need_new_partners2.print_agents()

                if len(partner._partners) >= partner._num_sex_partners :
                    #print "Partners", partner._partners
                    need_new_partners2.remove_agent(partner)

                if len(agent._partners) >= agent._num_sex_partners :
                    need_new_partners2.remove_agent(agent)
                    break


                missed_counter = 0
                partner = None

            else:
                #print "Missed pass attempt",missed_counter
                missed_counter += 1

            if missed_counter > 1:
                #print "No partner matches found"
                need_new_partners2.remove_agent(agent)
                noMatch += 1
                missed_counter = 0
                break

    # The remaining partnerless people can remain partnerless :)
    print "\n\t\t-COULDNT MATCH",noMatch,"AGENTS IN NEED \t---"
    #self.totalAgentClass.print_agents()




def reset_partner_count(self):
    """
    Reset partner count for method assess_interaction_distribution
    """

    # set ND partner count to zero for the next time step
    self.tmp_ND_NumPartners_Count = {}
    self.tmp_NIDU_NumPartners_Count = {}
    self.tmp_IDU_NumPartners_Count = {}
    self.tmp_MSM_NumPartners_Count = {}
