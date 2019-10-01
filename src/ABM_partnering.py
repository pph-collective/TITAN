#!/usr/bin/env python
# encoding: utf-8

# Imports
import random
from copy import deepcopy, copy
import os
import time

from functools import wraps
import numpy as np

try:
    from .HIVABM_Population import PopulationClass, print_population
except ImportError:
    raise ImportError("Can't import PopulationClass")

# IS THIS USED? CIRCULAR REFERENCES #REVIEW - delete
# try:
#     from .ABM_core import *
# except ImportError:
#     raise ImportError("Can't import PopulationClass")

try:
    from .agent import *
except ImportError:
    raise ImportError("Can't import Agent class")

from . import params


def update_partner_assignments(self, partnerTurnover, graph, agent=None):
    # Now create partnerships until available partnerships are out
    if agent:
        partner = get_partner(self, agent, self.All_agentSet)

        if partner:
            duration = get_partnership_duration(self, agent)
            tmp_relationship = Relationship(agent, partner, "MSM", "SE", duration)
            agent.bond(partner, tmp_relationship)
            self.Relationships.add_agent(tmp_relationship)
            graph.add_edge(tmp_relationship._ID1, tmp_relationship._ID2)
    else:
        EligibleAgents = self.All_agentSet
        noMatch = 0
        for agent in EligibleAgents.iter_agents():
            acquirePartnerProb = (
                params.cal_SexualPartScaling * partnerTurnover * (agent._mean_num_partners / (12.0))
            )
            if np.random.uniform(0, 1) < acquirePartnerProb:
                partner = get_partner(self, agent, self.All_agentSet)

                if partner:
                    duration = get_partnership_duration(self, agent)
                    tmp_relationship = Relationship(agent, partner, "MSM", "SE", duration)

                    agent.bond(partner, tmp_relationship)
                    self.Relationships.add_agent(tmp_relationship)
                    graph.add_edge(tmp_relationship._ID1, tmp_relationship._ID2)
                else:
                    noMatch += 1
                    graph.add_node(agent)
            else:
                graph.add_node(agent)

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
    shortlist_NNP = need_new_partners
    agent_race_type = agent._race
    agent_sex_type = agent._SO
    agent_drug_type = agent._DU
    RandomPartner = None

    if agent_drug_type == "IDU":
        if random.random() < 0.8:
            # choose from IDU agents
            try:
                RandomPartner = get_random_IDU_partner(self, agent, shortlist_NNP)
            except:
                print("No IDU matches")
                get_random_sex_partner(self, agent, shortlist_NNP)
        else:
            get_random_sex_partner(self, agent, shortlist_NNP)
    elif agent_drug_type in ("NDU", "NIDU"):
        if params.flag_AssortativeMix:
            if random.random() < params.DemographicParams[agent_race_type]["ALL"]["AssortMixCoeff"]:
                RandomPartner = get_assort_sex_partner(self, agent, shortlist_NNP)
                if not RandomPartner and params.AssortMixCoeff <= 1.0:
                    RandomPartner = get_random_sex_partner(self, agent, shortlist_NNP)
            else:
                RandomPartner = get_random_sex_partner(self, agent, shortlist_NNP)
        else:
            RandomPartner = get_random_sex_partner(self, agent, shortlist_NNP)
    else:
        raise ValueError("Check method _get_partners(). Agent not caught!")

    if RandomPartner == agent:
        return None
    else:
        return RandomPartner


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

    # todo: Make the random agent never return the agent or any of their partners
    if agent_drug_type not in ["IDU"]:
        raise ValueError("Invalid drug type! %s" % str(agent_drug_type))
    else:
        RandomPartner = random.choice([ptn for ptn in need_new_partners._subset["DU"]._subset["IDU"]._members if ptn not in agent._partners])
        if RandomPartner in agent._partners or RandomPartner == agent:
            RandomPartner = None

    if RandomPartner:
        return RandomPartner
    else:
        return None

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

    def getPartnerBin(agent):
        testRand = random.random()
        i = 1
        pMatch = params.mixingMatrix[agent._ageBin][i]

        if params.flag_AgeAssortMix:
            while True:
                if testRand <= pMatch:
                    return i
                else:
                    i += 1
                    pMatch += params.mixingMatrix[agent._ageBin][i]
                if i == 5:
                    return i
        else:
            i = random.randrange(1, 6)
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
    # todo: Make the random agent never return the agent or any of their partners
    assert agent_sex_type in ["HM", "HF", "MSM", "WSW", "MTF"]

    eligPartnerType = params.DemographicParams[agent_race_type][agent_sex_type][
        "EligSE_PartnerType"
    ][0]

    if params.AssortMixType == "Age":
        randomK_sample = random.sample(
            need_new_partners._subset["MSM"]._members, params.cal_ptnrSampleDepth
        )
        ageBinPick = getPartnerBin(agent)
        while True:
            availableParts = [ag for ag in randomK_sample if ag not in agent._partners]
            RandomPartner = random.choice([ag for ag in availableParts if ag._ageBin == ageBinPick])
            break

    # else if picking using race mix
    elif params.AssortMixType == "Race":
        samplePop = [
            tmpA
            for tmpA in need_new_partners._subset["SO"]._subset[eligPartnerType]._members
            if (tmpA._race == agent._race and tmpA not in agent._partners)
        ]
        try:
            randomK_sample = random.sample(samplePop, params.cal_ptnrSampleDepth)
        except:
            randomK_sample = samplePop
        while True:
            RandomPartner = random.choice(samplePop)
            break

    elif params.AssortMixType == "Client":
        if agent._race == "WHITE":
            samplePop = [
                tmpA
                for tmpA in need_new_partners._subset["SO"]._subset[eligPartnerType]._members
                if (tmpA._race == "WHITE" and tmpA not in agent._partners)
            ]
            try:
                randomK_sample = random.sample(samplePop, params.cal_ptnrSampleDepth)
            except:
                randomK_sample = samplePop
        else:
            samplePop = [
                tmpA
                for tmpA in need_new_partners._subset["SO"]._subset[eligPartnerType]._members
                if (tmpA._race == "WHITE" and tmpA._everhighrisk_bool and tmpA not in agent._partners)
            ]
            try:
                randomK_sample = random.sample(samplePop, params.cal_ptnrSampleDepth)
            except:
                randomK_sample = samplePop
        while True:
            RandomPartner = random.choice(samplePop)
            break

    elif params.AssortMixType == "HR":
        samplePop = [
            tmpA
            for tmpA in need_new_partners._subset["SO"]._subset[eligPartnerType]._members
            if (tmpA._everhighrisk_bool and tmpA not in agent._partners)
        ]
        if samplePop:
            try:
                randomK_sample = random.sample(samplePop, params.cal_ptnrSampleDepth)
            except:
                randomK_sample = samplePop

            while True:
                RandomPartner = random.choice(samplePop)
                break

    if RandomPartner == None or RandomPartner in agent._partners or RandomPartner == agent:
        RandomPartner = None
    else:
        try:
            RandomPartner = random.choice(tempList)
        except:
            pass

    if RandomPartner:
        return RandomPartner
    else:
        pass


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

    eligPtnType = params.DemographicParams[agent_race_type][agent_sex_type]["EligSE_PartnerType"][0]

    partnerPool2 = need_new_partners._subset["SO"]._subset[eligPtnType]._members

    RandomPartner = random.choice(partnerPool2)

    if agent_sex_type not in params.agentSexTypes:
        raise ValueError("Invalid sex type! %s" % str(agent_sex_type))
    else:
        pass

    if RandomPartner:
        if (RandomPartner in agent._partners) or (RandomPartner == agent):
            pass
        else:
            assert sex_possible(
                self, agent._SO, RandomPartner._SO
            ), "Sex no possible between agents! ERROR 441"
            return RandomPartner
    else:
        pass

#REVIEW is this redundant with ABM_Core._sex_possible ? - replace this with _sex_possible, delete _sex_possible and update references
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
    if agent_sex_type not in ["HM", "HF", "MSM", "WSW", "MTF", "MSW"]:
        raise ValueError("Invalid agent_sex_type! %s" % str(agent_sex_type))
    if partner_sex_type not in ["HM", "HF", "MSM", "WSW", "MTF", "MSW"]:
        raise ValueError("Invalid partner_sex_type! %s" % str(partner_sex_type))

    # Sex possible
    if agent_sex_type == "HM" and partner_sex_type in ["HF", "WSW", "MTF"]:
        SexPossible = True
    elif agent_sex_type == "MSM" and partner_sex_type in ["MSM", "WSW", "HF", "MTF", "MSW"]:
        SexPossible = True
    elif agent_sex_type == "WSW" and partner_sex_type in ["MSM", "WSW", "HM"]:
        SexPossible = True
    elif agent_sex_type == "HF" and partner_sex_type in ["HM", "MSM"]:
        SexPossible = True
    elif agent_sex_type == "MTF" and partner_sex_type in ["HM", "MSM"]:
        SexPossible = True
    elif agent_sex_type == "MSW" and partner_sex_type in ["MSM", "WSW", "HF", "MTF", "MSW"]:
        SexPossible = True
    else:
        SexPossible = False

    if agent_sex_type == "HM" and partner_sex_type == "HM" and SexPossible:
        raise ValueError("Check _sex_possible method!")

    return SexPossible


def get_partnership_duration(self, agent):
    """
    :Purpose:
        Get number of partners for a agent.
        Drawn from Poisson distribution.

    :Input:
        agent : Agent

    :Output:
        NumPartners : int
        Zero partners possible.
    """
    # Check input
    agent_drug_type = agent._DU
    agent_sex_type = agent._SO
    agent_race_type = agent._race

    # Drug type
    if agent_drug_type not in ["IDU", "NIDU", "NDU"]:
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
    if agent_race_type == "BLACK" and params.model == "MSW":
        MSWsexualDurations = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
        MSWsexualDurations[1] = {"p_value": (0.27 + 0.22), "min": 1, "max": 6}
        MSWsexualDurations[2] = {"p_value": (0.09 + 0.262 + 0.116), "min": 7, "max": 12}
        MSWsexualDurations[3] = {"p_value": (0.09 + 0.09), "min": 13, "max": 24}
        MSWsexualDurations[4] = {"p_value": (0.09 + 0.09 + 0.07), "min": 25, "max": 36}
        MSWsexualDurations[5] = {"min": 37, "max": 48}

        if diceroll < MSWsexualDurations[1]["p_value"]:
            dur_bin = 1
        elif diceroll < MSWsexualDurations[2]["p_value"]:
            dur_bin = 2
        elif diceroll < MSWsexualDurations[3]["p_value"]:
            dur_bin = 3
        elif diceroll < MSWsexualDurations[4]["p_value"]:
            dur_bin = 4
        else:
            dur_bin = 5

        duration = random.randrange(
            MSWsexualDurations[dur_bin]["min"], MSWsexualDurations[dur_bin]["max"], 1
        )

    else:
        if diceroll < params.sexualDurations[1]["p_value"]:
            dur_bin = 1
        elif diceroll < params.sexualDurations[2]["p_value"]:
            dur_bin = 2
        elif diceroll < params.sexualDurations[3]["p_value"]:
            dur_bin = 3
        elif diceroll < params.sexualDurations[4]["p_value"]:
            dur_bin = 4
        else:
            dur_bin = 5

        duration = random.randrange(
            params.sexualDurations[dur_bin]["min"], params.sexualDurations[dur_bin]["max"], 1
        )

    return duration
