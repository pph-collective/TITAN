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
#import random
from copy import deepcopy, copy
import os
import time
#import PyQt4

#import scipy.sparse as spsp
# from scipy.stats import binom
# from scipy.stats import poisson
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
    # print "Update partner random start ",random.randint(0,1000)
    if agent:
        partner = get_partner(self, agent, self.All_agentSet)

        if partner:
            #print "Agent %d found partner %d!"%(agent.get_ID(), partner.get_ID())
            duration = get_partnership_duration(self, agent)
            tmp_relationship = Relationship(agent, partner, "MSM", "SE", duration)
            agent.bond(partner, tmp_relationship)
            self.Relationships.add_agent(tmp_relationship)
            self.G.add_edge(tmp_relationship._ID1, tmp_relationship._ID2)
    else:
        EligibleAgents = self.All_agentSet#._subset["HIV"].iter_agents()
        noMatch = 0
        for agent in EligibleAgents.iter_agents():
            #print len(agent._partners)
            acquirePartnerProb = params.cal_SexualPartScaling*partnerTurnover*(agent._mean_num_partners / (12.0))
            #if agent._highrisk_bool:print acquirePartnerProb
            if np.random.uniform(0, 1) < acquirePartnerProb:
                partner = get_partner(self, agent, self.All_agentSet)

                if partner:
                    #print "Agent %d found partner %d!"%(agent.get_ID(), partner.get_ID())

                    duration = get_partnership_duration(self, agent)
                    tmp_relationship = Relationship(agent, partner, "MSM", "SE", duration)

                    agent.bond(partner, tmp_relationship)
                    self.Relationships.add_agent(tmp_relationship)
                    self.G.add_edge(tmp_relationship._ID1, tmp_relationship._ID2)

                    # ADD RELATIONSHIP EDGE TO GRAPH G of NetworkGraph
                    #print "%d/%d partnets found for agent %d"%(len(agent._partners), agent._num_sex_partners, agent.get_ID())
                    #print "%d/%d partnets found for partner %d"%(len(partner._partners), partner._num_sex_partners, partner.get_ID())

                else:
                    # print "Missed pass attempt",noMatch
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
        RandNumCont = RandNumCont * params.cal_SexualPartScaling


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
    agent_race_type = agent._race
    agent_sex_type = agent._SO
    agent_drug_type = agent._DU
    RandomPartner = None


    #print("Finding partner for agent", agent._ID, agent_sex_type, agent_drug_type)
    if agent_drug_type == 'IDU':
        if random.random() < 0.8:
            # choose from IDU agents
            try: RandomPartner = get_random_IDU_partner(self, agent, shortlist_NNP)
            except:
                print("No IDU matches")
                get_random_sex_partner(self, agent, shortlist_NNP)
            #print "\tReturned: %s" % RandomPartner
        else:
            get_random_sex_partner(self, agent, shortlist_NNP)
    elif agent_drug_type in ('NDU','NIDU'):
        if params.flag_AssortativeMix:
            if random.random() < params.DemographicParams[agent_race_type]['ALL']['AssortMixCoeff']:
                RandomPartner = get_assort_sex_partner(self, agent, shortlist_NNP)
                if not RandomPartner and params.AssortMixCoeff <= 1.0:
                    RandomPartner = get_random_sex_partner(self, agent, shortlist_NNP)
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
        RandomPartner = random.choice(need_new_partners._subset["DU"]._subset["IDU"]._members)
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
            Assortive_intersection = list(set(self.Race_WHITE_agentSet).intersection(intersection))
            if Assortive_intersection == []: print("Couldnt assortive mix (W), picking suitable agent")
            else: return random.choice(Assortive_intersection)
        elif agent_race_type == 'BLACK':
            Assortive_intersection = list(set(self.Race_BLACK_agentSet).intersection(intersection))
            if Assortive_intersection == []:
                print("Couldnt assortive mix (B), picking suitable agent")
            else:
                #print Assortive_intersection
                return random.choice(Assortive_intersection)
        if intersection == []: return None
        else: print("NO PATNAS")#return random.choice(intersection)

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

    eligPartnerType = params.DemographicParams[agent_race_type][agent_sex_type]['EligSE_PartnerType'][0]

    if params.AssortMixType == 'Age':
        randomK_sample = random.sample(need_new_partners._subset["MSM"]._members, params.cal_ptnrSampleDepth)
        ageBinPick = getPartnerBin(agent)
        while True:
            RandomPartner = random.choice([ag for ag in randomK_sample if ag._ageBin == ageBinPick])
            break

    #else if picking using race mix
    elif params.AssortMixType == 'Race':
        samplePop = [tmpA for tmpA in need_new_partners._subset['SO']._subset[eligPartnerType]._members if tmpA._race == agent._race]
        try:
            randomK_sample = random.sample(samplePop ,params.cal_ptnrSampleDepth)
        except:
            randomK_sample = samplePop
        while True:
            RandomPartner = random.choice(samplePop)
            break

    elif params.AssortMixType == 'Client':
        if agent._race == 'WHITE':
            samplePop = [tmpA for tmpA in need_new_partners._subset['SO']._subset[eligPartnerType]._members if (tmpA._race == 'WHITE')]
            try:
                randomK_sample = random.sample(samplePop ,params.cal_ptnrSampleDepth)
            except:
                randomK_sample = samplePop
        else:
            samplePop = [tmpA for tmpA in need_new_partners._subset['SO']._subset[eligPartnerType]._members if (tmpA._race == 'WHITE' and tmpA._everhighrisk_bool)]
            try:
                randomK_sample = random.sample(samplePop ,params.cal_ptnrSampleDepth)
            except:
                randomK_sample = samplePop
        while True:
            RandomPartner = random.choice(samplePop)
            break

    elif params.AssortMixType == 'HR':
        samplePop = [tmpA for tmpA in need_new_partners._subset['SO']._subset[eligPartnerType]._members if tmpA._everhighrisk_bool]
        if samplePop:
            try:
                randomK_sample = random.sample(samplePop ,params.cal_ptnrSampleDepth)
            except:
                randomK_sample = samplePop

            while True:
                RandomPartner = random.choice(samplePop)
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

    agent_sex_type = agent._SO
    agent_race_type = agent._race
    agent_drug_type = agent._DU

    RandomPartner = None

    # partnerPool = set()
    # AssortMix = False
    # for eligPtnType in params.DemographicParams[agent_race_type][agent_sex_type]['EligSE_PartnerType']:
    #     partnerPool = (need_new_partners._subset["SO"]._subset[eligPtnType]._members)
    # print type(partnerPool)
    # print type(need_new_partners)
    # print type(need_new_partners._subset[eligPtnType]._members)
    eligPtnType = params.DemographicParams[agent_race_type][agent_sex_type]['EligSE_PartnerType'][0]
    #partnerPool = list(need_new_partners._subset["SO"]._subset[eligPtnType]._members)
    #partnerPool = list(need_new_partners._subset["SO"]._subset[eligPtnType]._members- set(agent._partners) - set([agent]))
    #removeFromPool = set(agent._partners) - set([agent])
    partnerPool2 = need_new_partners._subset["SO"]._subset[eligPtnType]._members
    # print(partnerPool2)
    # RandomPartner = (need_new_partners._subset["SO"]._subset[eligPtnType]._members- set(agent._partners) - set([agent])))
    RandomPartner = random.choice(partnerPool2)
    # print type(partnerPool)
    # print type(RandomPartner)

    if agent_sex_type not in params.agentSexTypes:
        raise ValueError("Invalid sex type! %s"%str(agent_sex_type))
    else:
        pass
        # RandomPartner = random.choice(tuple(partnerPool))
        # RandomPartner = random.sample((partnerPool), 1)
        #RandomPartner = partnerPool.pop()
        #partnerPool.add(RandomPartner)


    if RandomPartner:
        if (RandomPartner in agent._partners) or (RandomPartner == agent):
            pass
        # assert(sex_possible(self,agent._SO, RandomPartner[0]._SO)),"Sex no possible between agents! ERROR 441"
        # return RandomPartner[0]
        else:
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
    if agent_sex_type not in ['HM', 'HF', 'MSM', 'WSW','MTF','MSW']:
        raise ValueError("Invalid agent_sex_type! %s" % str(agent_sex_type))
    if partner_sex_type not in ['HM', 'HF', 'MSM', 'WSW','MTF','MSW']:
        raise ValueError("Invalid partner_sex_type! %s" % str(
            partner_sex_type))

    # Sex possible
    if agent_sex_type == 'HM' and partner_sex_type in ['HF', 'WSW', 'MTF']:
        SexPossible = True
    #elif partner_sex_type == 'HM' and agent_sex_type in ['HF', 'WSW']:
    #    SexPossible = True
    elif agent_sex_type == 'MSM' and partner_sex_type in ['MSM', 'WSW', 'HF', 'MTF', 'MSW']:
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
    elif agent_sex_type == 'MSW' and partner_sex_type in ['MSM', 'WSW', 'HF', 'MTF', 'MSW']:
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
    agent_race_type = agent._race

    # Drug type
    if agent_drug_type not in ['IDU', 'NIDU', 'NDU']:
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
    if agent_race_type == "BLACK" and params.model == 'MSW':
        MSWsexualDurations = {1:{}, 2:{}, 3:{}, 4:{}, 5:{}}
        MSWsexualDurations[1] = {'p_value':(0.27 + 0.22), 'min':1, 'max':6}
        MSWsexualDurations[2] = {'p_value':(0.09 + 0.262 + 0.116), 'min':7, 'max':12}
        MSWsexualDurations[3] = {'p_value':(0.09 + 0.09), 'min':13, 'max':24}
        MSWsexualDurations[4] = {'p_value':(0.09 + 0.09 + 0.07), 'min':25, 'max':36}
        MSWsexualDurations[5] = {'min':37, 'max':48}

        if diceroll < MSWsexualDurations[1]['p_value']:
            dur_bin = 1
        elif diceroll < MSWsexualDurations[2]['p_value']:
            dur_bin = 2
        elif diceroll < MSWsexualDurations[3]['p_value']:
            dur_bin = 3
        elif diceroll < MSWsexualDurations[4]['p_value']:
            dur_bin = 4
        else:
            dur_bin = 5

        duration = random.randrange(MSWsexualDurations[dur_bin]['min'], MSWsexualDurations[dur_bin]['max'], 1)

    else:
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



def reset_partner_count(self):
    """
    Reset partner count for method assess_interaction_distribution
    """

    # set ND partner count to zero for the next time step
    self.tmp_ND_NumPartners_Count = {}
    self.tmp_NIDU_NumPartners_Count = {}
    self.tmp_IDU_NumPartners_Count = {}
    self.tmp_MSM_NumPartners_Count = {}
