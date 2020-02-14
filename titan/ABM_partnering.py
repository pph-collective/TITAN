#!/usr/bin/env python
# encoding: utf-8

# Imports
import random
from typing import Sequence, List, Dict, Optional, TypeVar

from dotmap import DotMap  # type: ignore

from . import probabilities as prob
from .agent import Agent, AgentSet


def get_partner(
    agent: Agent, all_agent_set: AgentSet, params: DotMap
) -> Optional[Agent]:
    """
    :Purpose:
        Get partner for agent.

    :Input:
        agent : Agent

        all_agent_set: AgentSet of partners to pair with

    :Output:
        partner: new partner
    """
    agent_drug_type = agent.drug_use
    RandomPartner = None

    if agent_drug_type == "Inj":
        if random.random() < 0.8:
            # choose from PWID agents
            RandomPartner = get_random_PWID_partner(agent, all_agent_set)

        # either didn't try to get PWID partner, or failed to get PWID partner
        if RandomPartner is None:
            get_random_sex_partner(agent, all_agent_set, params)
    elif agent_drug_type in ("None", "NonInj"):
        if params.features.assort_mix and (
            random.random() < params.demographics[agent.race].assort_mix.coefficient
        ):
            RandomPartner = get_assort_sex_partner(agent, all_agent_set, params)

            # try again with random sex partner is assort mix percent not 100%
            if RandomPartner is None and params.assort_mix.coefficient <= 1.0:
                RandomPartner = get_random_sex_partner(agent, all_agent_set, params)
        else:
            RandomPartner = get_random_sex_partner(agent, all_agent_set, params)
    else:
        raise ValueError("Check method _get_partners(). Agent not caught!")

    if RandomPartner == agent:
        return None
    else:
        return RandomPartner


def get_random_PWID_partner(agent: Agent, all_agent_set: AgentSet) -> Optional[Agent]:
    """
    :Purpose:
        Get a random partner which is sex compatible

    :Input:
        agent: Agent
        all_agent_set: AgentSet of partners to pair with

    :Output:
        partner : Agent or None

    """
    agent_drug_type = agent.drug_use
    assert agent_drug_type == "Inj", "Agent's drug type must be Inj"

    RandomPartner = None

    RandomPartner = safe_random_choice(
        [
            ptn
            for ptn in all_agent_set.subset["DU"].subset["Inj"].members
            if ptn not in agent.partners and ptn != agent
        ]
    )

    return RandomPartner


def get_assort_sex_partner(
    agent: Agent, all_agent_set: AgentSet, params: DotMap
) -> Optional[Agent]:
    """
    :Purpose:
        Get a random partner which is sex compatible and fits assortativity constraints

    :Input:
        agent: int
        all_agent_set: AgentSet of partners to pair with

    :Output:
        partner : Agent or None

    """

    RandomPartner = None

    assert agent.so in params.classes.sex_types

    eligible_partners = [
        partner
        for partner in all_agent_set
        if sex_possible(agent.so, partner.so, params)
    ]

    if params.assort_mix.type == "Race":
        samplePop = [
            tmpA
            for tmpA in eligible_partners
            if (tmpA.race == agent.race and tmpA not in agent.partners)
        ]

    elif params.assort_mix.type == "Client":
        if agent.race == "WHITE":
            samplePop = [
                tmpA
                for tmpA in eligible_partners
                if (tmpA.race == "WHITE" and tmpA not in agent.partners)
            ]
        else:
            samplePop = [
                tmpA
                for tmpA in eligible_partners
                if (
                    tmpA.race == "WHITE"
                    and tmpA.high_risk_ever
                    and tmpA not in agent.partners
                )
            ]

    elif params.assort_mix.type == "high_risk":
        samplePop = [
            tmpA
            for tmpA in eligible_partners
            if (tmpA.high_risk_ever and tmpA not in agent.partners)
        ]

    RandomPartner = safe_random_choice(samplePop)

    # partner can't be existing parter or agent themself
    if RandomPartner in agent.partners or RandomPartner == agent:
        RandomPartner = None

    return RandomPartner


def get_random_sex_partner(
    agent: Agent, all_agent_set: AgentSet, params: DotMap
) -> Optional[Agent]:
    """
    :Purpose:
        Get a random partner which is sex compatible

    :Input:
        agent: Agent
        all_agent_set: list of available partners (AgentSet)

    :Output:
        partner : Agent or None

    """
    random_partner = None

    elig_partner_pool = [
        partner
        for partner in all_agent_set.members
        if sex_possible(agent.so, partner.so, params)
    ]

    random_partner = safe_random_choice(elig_partner_pool)

    if (random_partner in agent.partners) or (random_partner == agent):
        random_partner = None

    return random_partner


def sex_possible(agent_sex_type: str, partner_sex_type: str, params: DotMap) -> bool:
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
        "HM": ["HF", "MTF"],
        "MSM": ["MSM", "MTF"],
        "WSW": ["WSW", "MTF"],
        "HF": ["HM"],
        "MTF": ["WSW", "HM", "MSM"],
    }

    # Check input
    if agent_sex_type not in params.classes.sex_types:
        raise ValueError("Invalid agent_sex_type! %s" % str(agent_sex_type))
    if partner_sex_type not in params.classes.sex_types:
        raise ValueError("Invalid partner_sex_type! %s" % str(partner_sex_type))

    return (agent_sex_type in st[partner_sex_type]) and (
        partner_sex_type in st[agent_sex_type]
    )


def get_partnership_duration(agent: Agent, params: DotMap) -> int:
    """
    :Purpose:
        Get duration of a relationship
        Drawn from Poisson distribution.

    :Input:
        agent : Agent

    :Output:
        NumPartners : int
        Zero partners possible.
    """
    # Check input
    agent_drug_type = agent.drug_use
    agent_sex_type = agent.so
    agent_race_type = agent.race

    # Drug type
    if agent_drug_type not in ["Inj", "NonInj", "None"]:
        raise ValueError("Invalid drug type! %s" % str(agent_drug_type))
    # Sex type
    if agent_sex_type not in params.classes.sex_types:
        raise ValueError("Invalid sex type! %s" % str(agent_sex_type))

    diceroll = random.random()

    # Length of relationship (months)a
    # <1 1,679 32.3% 566 17.7 1,113 55.8
    # 1–6 1,359 26.2% 929 29.0 430 21.6
    # 7–12 604 11.6% 459 14.4 145 7.3
    # 13–24 628 12.1% 480 15.0 148 7.4
    # 25–36 309 6.0% 264 8.3 45 2.3
    # >37 614 11.8% 501 15.7 113 5.7
    if (
        agent_race_type == "BLACK" and params.model.features.msmw
    ):  # TO_REVIEW is this right?
        dur_bin = 5
        for i in range(1, 5):
            if diceroll < prob.MSWsexualDurations[i]["p_value"]:
                dur_bin = i
                break

        duration = random.randrange(
            prob.MSWsexualDurations[dur_bin]["min"],
            prob.MSWsexualDurations[dur_bin]["max"],
            1,
        )

    else:
        dur_bin = 5
        for i in range(1, 5):
            if diceroll < params.partnership.sex.duration[i].prob:
                dur_bin = i
                break

        duration = random.randrange(
            params.partnership.sex.duration[dur_bin].min,
            params.partnership.sex.duration[dur_bin].max,
            1,
        )

    return duration


T = TypeVar("T")


def safe_random_choice(seq: Sequence[T]) -> Optional[T]:
    """
    Return None or a random choice
    """
    if seq:
        return random.choice(seq)
    else:
        return None
