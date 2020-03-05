#!/usr/bin/env python
# encoding: utf-8

# Imports
import random
from typing import Optional, Tuple, Set

from dotmap import DotMap  # type: ignore

from .agent import Agent, AgentSet
from . import utils


def select_partner(
    agent: Agent, need_partners: Set, params: DotMap, rand_gen
) -> Tuple[Optional[Agent], str]:
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
    partner = None
    partner_set: Set[Agent] = set(need_partners)
    eligible_partner_set = partner_set - set(agent.partners) - {agent}
    RandomPartner: Optional[Agent]

    def bondtype(bond_dict):
        bonds = {"type": [], "prob": []}
        for value in bond_dict[agent.race]:
            bonds["type"].append(value["type"])
            bonds["prob"].append(value["probability"])
        print(bonds)
        bonded_type = rand_gen.choices(bonds["type"], weights=bonds["prob"], k=1)
        return bonded_type

    def assort(eligible_partner_list, assort_params):
        if rand_gen.random() < assort_params["probability"]:
            eligible_partners = {
                partner
                for partner in eligible_partner_list
                if (
                    getattr(partner, assort_params["assort_type"])
                    == assort_params["partner_type"]
                )
            }
        else:
            eligible_partners = {
                partner
                for partner in eligible_partner_list
                if (
                    getattr(partner, assort_params["assort_type"])
                    != assort_params["partner_type"]
                )
            }
        return eligible_partners

    if params.features.assort:
        for assort_types in params.assort_mix.assortativity.assort_type:
            if getattr(agent, assort_types.assort_type) == assort_types["agent_type"]:
                eligible_partner_set = assort(eligible_partner_set, assort_types)

    if agent_drug_type == "Inj":
        agent_bond = bondtype(params.partnering.bond)
    else:
        agent_bond = bondtype(params.partnering.bond)

    if "injection" in agent_bond:
        eligible_partner_set = {
            partner for partner in eligible_partner_set if partner.drug_use == "Inj"
        }
    if "sexual" in agent_bond:
        eligible_partner_set = {
            partner
            for partner in eligible_partner_set
            if sex_possible(agent.so, partner.so, params)
        }
    if "social" in agent_bond:
        eligible_partner_set = eligible_partner_set

    if eligible_partner_set:
        random_partner = rand_gen.choices(eligible_partner_set, k=1)[0]
    else:
        random_partner = None

    return random_partner, agent_bond


def get_random_pwid_partner(
    agent: Agent, all_agent_set: AgentSet, rand_gen
) -> Optional[Agent]:
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

    partner = None

    partner = utils.safe_random_choice(
        [
            p
            for p in all_agent_set.subset["DU"].subset["Inj"].members
            if p not in agent.partners and p != agent
        ],
        rand_gen,
    )

    return partner


def get_assort_sex_partner(
    agent: Agent, all_agent_set: AgentSet, params: DotMap, rand_gen
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

    partner = None

    assert agent.so in params.classes.sex_types

    eligible_partners = [
        p for p in all_agent_set if sex_possible(agent.so, p.so, params)
    ]
    return partner


def get_random_sex_partner(
    agent: Agent, all_agent_set: AgentSet, params: DotMap, rand_gen
) -> Optional[Agent]:
    """
    :Purpose:
        Get a random partner which is sex compatible

    :Input:
        agent: Agent
        all_agent_set: list of available partners (AgentSet)
        params: model parameters

    :Output:
        partner : Agent or None

    """
    partner = None

    eligible_partners = [
        p for p in all_agent_set.members if sex_possible(agent.so, p.so, params)
    ]

    partner = utils.safe_random_choice(eligible_partners, rand_gen)

    if (partner in agent.partners) or (partner == agent):
        partner = None

    return partner


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

    # Check input
    if agent_sex_type not in params.classes.sex_types:
        raise ValueError("Invalid agent_sex_type! %s" % str(agent_sex_type))
    if partner_sex_type not in params.classes.sex_types:
        raise ValueError("Invalid partner_sex_type! %s" % str(partner_sex_type))

    agent_match = (
        agent_sex_type in params.classes.sex_types[partner_sex_type].sleeps_with
    )
    partner_match = (
        partner_sex_type in params.classes.sex_types[agent_sex_type].sleeps_with
    )

    return agent_match and partner_match


def get_partnership_duration(agent: Agent, params: DotMap, rand_gen) -> int:
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

    # Length of relationship (months)a
    # <1 1,679 32.3% 566 17.7 1,113 55.8
    # 1–6 1,359 26.2% 929 29.0 430 21.6
    # 7–12 604 11.6% 459 14.4 145 7.3
    # 13–24 628 12.1% 480 15.0 148 7.4
    # 25–36 309 6.0% 264 8.3 45 2.3
    # >37 614 11.8% 501 15.7 113 5.7
    diceroll = rand_gen.random()
    dur_bin = 5
    for i in range(1, 5):
        if diceroll < params.partnership.sex.duration[i].prob:
            dur_bin = i
            break

    duration = rand_gen.randrange(
        params.partnership.sex.duration[dur_bin].min,
        params.partnership.sex.duration[dur_bin].max,
        1,
    )

    return duration
