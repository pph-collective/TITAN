#!/usr/bin/env python
# encoding: utf-8

# Imports
import random
from typing import Sequence, Optional, TypeVar

from . import params  # type: ignore
from . import probabilities as prob
from .agent import Agent, Agent_set


def get_partner(agent: Agent, all_agent_set: Agent_set) -> tuple:
    """
    :Purpose:
        Get partner for agent.

    :Input:
        agent : Agent

        all_agent_set: AgentSet of partners to pair with

    :Output:
        partner: new partner
    """
    agent_drug_type = agent._DU
    eligible_partners = all_agent_set._members
    # eligible_partners.remove_agent(agent)

    def bondtype(bond_dict):
        pvalue = 0.0
        bond_probability = random.random()
        bonded_type = "sexualOnly"
        for key, values in bond_dict.items():
            pvalue += values["probability"]
            if bond_probability < pvalue:
                bonded_type = values["type"]
                break
        return bonded_type

    def assort(eligible_partners):
        for assort_key, assort_values in params.assortative_mixing.items():
            if getattr(agent, assort_values["type"]) == assort_values["agent_type"]:
                assert (
                    assort_values["probability"] != 0.0
                ), "Cannot assort with probability of 0!"
                if (
                    random.random() < assort_values["probability"]
                ):  # if roll to assortatively mix, only that class is
                    # your eligible partners
                    eligible_partners = [
                        tmpA
                        for tmpA in eligible_partners
                        if (
                            getattr(tmpA, assort_values["type"])
                            == assort_values["type"]
                        )
                    ]
                else:  # if you aren't chosen to assortative mix,
                    # you should be choosing people who don't meet that criteria
                    eligible_partners = [
                        tmpA
                        for tmpA in eligible_partners
                        if (
                            getattr(tmpA, assort_values["type"])
                            != assort_values["type"]
                        )
                    ]
        return eligible_partners

    if (
        agent_drug_type == "IDU"
    ):  # agent bond types are drug-use specific. Get the correct dict.
        agent_bond = bondtype(params.bond_type_probs_IDU)
    else:
        agent_bond = bondtype(params.bond_type_probs)

    if "injection" in agent_bond:
        eligible_partners = [
            ptn
            for ptn in all_agent_set._subset["DU"]._subset["IDU"]._members
            if ptn not in agent._partners and ptn != agent
        ]
    if "sexual" in agent_bond:
        eligible_partners = [
            ptn
            for ptn in all_agent_set.iter_agents()
            if sex_possible(agent._SO, ptn._SO)
            and ptn not in agent._partners
            and ptn != agent
        ]
    elif "social" in agent_bond:
        eligible_partners = [
            ptn
            for ptn in all_agent_set.iter_agents()
            if ptn not in agent._partners and ptn != agent
        ]

    assort(eligible_partners)
    RandomPartner = safe_random_choice(eligible_partners)
    return RandomPartner, agent_bond


def sex_possible(agent_sex_type: str, partner_sex_type: str) -> bool:
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


def get_partnership_duration(agent: Agent) -> int:
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


T = TypeVar("T")


def safe_random_choice(seq: Sequence[T]) -> Optional[T]:
    """
    Return None or a random choice
    """
    if seq:
        return random.choice(seq)
    else:
        return None
