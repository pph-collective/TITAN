#!/usr/bin/env python
# encoding: utf-8

# Imports
import random

from . import params
from . import probabilities as prob


def get_partner(agent, need_new_partners):
    """
    :Purpose:
        Get partner for agent.

    :Input:
        agent : Agent

        need_new_partners: AgentSet of partners to pair with

    :Output:
        partner: new partner
    """
    agent_race_type = agent._race
    agent_drug_type = agent._DU
    RandomPartner = None

    if agent_drug_type == "IDU":
        if random.random() < 0.8:
            # choose from IDU agents
            try:  # REVIEW perhaps get_random_IDU_partner should return none instead of throw exception
                RandomPartner = get_random_IDU_partner(agent, need_new_partners)
            except:
                print("No IDU matches")
                get_random_sex_partner(agent, need_new_partners)
        else:
            get_random_sex_partner(agent, need_new_partners)
    elif agent_drug_type in ("NDU", "NIDU"):
        if params.flag_AssortativeMix:
            if (
                random.random()
                < params.DemographicParams[agent_race_type]["ALL"]["AssortMixCoeff"]
            ):
                RandomPartner = get_assort_sex_partner(agent, need_new_partners)
                if not RandomPartner and params.AssortMixCoeff <= 1.0:
                    RandomPartner = get_random_sex_partner(agent, need_new_partners)
            else:
                RandomPartner = get_random_sex_partner(agent, need_new_partners)
        else:
            RandomPartner = get_random_sex_partner(agent, need_new_partners)
    else:
        raise ValueError("Check method _get_partners(). Agent not caught!")

    if RandomPartner == agent:
        return None
    else:
        return RandomPartner


def get_random_IDU_partner(agent, need_new_partners):
    """
    :Purpose:
        Get a random partner which is sex compatible

    :Input:
        agent: int
        need_new_partners: AgentSet of partners to pair with

    :Output:
        partner : int

    """
    agent_drug_type = agent._DU
    RandomPartner = None

    # REVIEW AssortMix never used
    AssortMix = False
    if random.random() < params.AssortMixCoeff:
        AssortMix = True

    # todo: Make the random agent never return the agent or any of their partners
    if agent_drug_type not in ["IDU"]:
        raise ValueError("Invalid drug type! %s" % str(agent_drug_type))
    else:
        RandomPartner = random.choice(
            [
                ptn
                for ptn in need_new_partners._subset["DU"]._subset["IDU"]._members
                if ptn not in agent._partners
            ]
        )
        if RandomPartner in agent._partners or RandomPartner == agent:
            RandomPartner = None

    return RandomPartner


def get_assort_sex_partner(agent, need_new_partners):
    """
    :Purpose:
        Get a random partner which is sex compatible and fits assortativity constraints

    :Input:
        agent: int
        need_new_partners: AgentSet of partners to pair with

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
    agent_race_type = agent._race

    RandomPartner = None
    tempList = []

    # REVIEW AssortMix never used - Sarah to reiew
    if random.random() < params.AssortMixCoeff:
        AssortMix = True
    else:
        AssortMix = False

    # todo: Make the random agent never return the agent or any of their partners
    assert agent_sex_type in [
        "HM",
        "HF",
        "MSM",
        "WSW",
        "MTF",
    ]  # REVIEW shouldn't this be the params. agent sex types? - switch this over

    # REVIEW in default params this is an empty list, in atlanta it's a len 1 list - why is this a list at all? (and is it eligible or ineligible)
    eligPartnerType = params.DemographicParams[agent_race_type][agent_sex_type][
        "EligSE_PartnerType"
    ][
        0
    ]  # Make this not a list - CHECK ALL OF THE SETTINGS

    if params.AssortMixType == "Age":
        randomK_sample = random.sample(
            need_new_partners._subset["MSM"]._members, params.cal_ptnrSampleDepth
        )
        ageBinPick = getPartnerBin(agent)
        while True:
            availableParts = [ag for ag in randomK_sample if ag not in agent._partners]
            RandomPartner = random.choice(
                [ag for ag in availableParts if ag._ageBin == ageBinPick]
            )
            break

    # else if picking using race mix
    elif params.AssortMixType == "Race":
        samplePop = [
            tmpA
            for tmpA in need_new_partners._subset["SO"]
            ._subset[eligPartnerType]
            ._members
            if (tmpA._race == agent._race and tmpA not in agent._partners)
        ]
        try:  # REVIEW refactor from try/except?
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
                for tmpA in need_new_partners._subset["SO"]
                ._subset[eligPartnerType]
                ._members
                if (tmpA._race == "WHITE" and tmpA not in agent._partners)
            ]
            try:  # REVIEW refactor from try/except?
                randomK_sample = random.sample(samplePop, params.cal_ptnrSampleDepth)
            except:
                randomK_sample = samplePop
        else:
            samplePop = [
                tmpA
                for tmpA in need_new_partners._subset["SO"]
                ._subset[eligPartnerType]
                ._members
                if (
                    tmpA._race == "WHITE"
                    and tmpA._everhighrisk_bool
                    and tmpA not in agent._partners
                )
            ]
            try:  # REVIEW refactor from try/except?
                randomK_sample = random.sample(samplePop, params.cal_ptnrSampleDepth)
            except:
                randomK_sample = samplePop
        while True:
            RandomPartner = random.choice(samplePop)
            break

    elif params.AssortMixType == "HR":
        samplePop = [
            tmpA
            for tmpA in need_new_partners._subset["SO"]
            ._subset[eligPartnerType]
            ._members
            if (tmpA._everhighrisk_bool and tmpA not in agent._partners)
        ]
        if samplePop:
            try:  # REVIEW refactor from try/except? (also generally repetitive code)
                randomK_sample = random.sample(samplePop, params.cal_ptnrSampleDepth)
            except:
                randomK_sample = samplePop

            while True:
                RandomPartner = random.choice(samplePop)
                break

    if (
        RandomPartner is None
        or RandomPartner in agent._partners
        or RandomPartner == agent
    ):
        RandomPartner = None
    else:
        try:
            RandomPartner = random.choice(tempList)
        except:  # REVIEW refactor from try/except?
            pass

    if RandomPartner:
        return RandomPartner
    else:
        pass


def get_random_sex_partner(agent, need_new_partners):
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

    RandomPartner = None

    eligPtnType = params.DemographicParams[agent_race_type][agent_sex_type][
        "EligSE_PartnerType"
    ][0]

    elig_partner_pool = need_new_partners._subset["SO"]._subset[eligPtnType]._members

    RandomPartner = random.choice(elig_partner_pool)

    if agent_sex_type not in params.agentSexTypes:
        raise ValueError("Invalid sex type! %s" % str(agent_sex_type))
    else:
        pass

    if RandomPartner:
        if (RandomPartner in agent._partners) or (RandomPartner == agent):
            pass
        else:
            assert sex_possible(
                agent._SO, RandomPartner._SO
            ), "Sex no possible between agents! ERROR 441"
            return RandomPartner
    else:
        pass


def sex_possible(agent_sex_type, partner_sex_type):
    """
    :Purpose:
    Determine if sex is possible.

    :Input:
    agent_sex_type : str

    partner_sex_type : str

    :Output:
    SexPossible : bool
    """

    # dictionary defining which sex types each sex type is compatible with
    st = {
        "HM": ["HF", "WSW", "MTF"],
        "MSM": ["MSM", "WSW", "HF", "MTF"],
        "WSW": ["MSM", "WSW", "HM"],
        "HF": ["HM", "MSM"],
        "MTF": ["HM", "MSM"],
    }

    # Check input
    if agent_sex_type not in params.agentSexTypes:
        raise ValueError("Invalid agent_sex_type! %s" % str(agent_sex_type))
    if partner_sex_type not in params.agentSexTypes:
        raise ValueError("Invalid partner_sex_type! %s" % str(partner_sex_type))

    return (agent_sex_type in st[partner_sex_type]) and (
        partner_sex_type in st[agent_sex_type]
    )


def get_partnership_duration(agent):
    """
    :Purpose:
        Get duration of a relationship #REVIEW sarah to look at why this is at the agent level
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

        if diceroll < prob.MSWsexualDurations[1]["p_value"]:
            dur_bin = 1
        elif diceroll < prob.MSWsexualDurations[2]["p_value"]:
            dur_bin = 2
        elif diceroll < prob.MSWsexualDurations[3]["p_value"]:
            dur_bin = 3
        elif diceroll < prob.MSWsexualDurations[4]["p_value"]:
            dur_bin = 4
        else:
            dur_bin = 5

        duration = random.randrange(
            prob.MSWsexualDurations[dur_bin]["min"],
            prob.MSWsexualDurations[dur_bin]["max"],
            1,
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
            params.sexualDurations[dur_bin]["min"],
            params.sexualDurations[dur_bin]["max"],
            1,
        )

    return duration
