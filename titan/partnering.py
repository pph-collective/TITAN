#!/usr/bin/env python
# encoding: utf-8

# Imports
from typing import Set, Optional, Tuple, Dict

from .agent import Agent, AgentSet
from . import utils
from .parse_params import ObjMap


def select_partner(
    agent: Agent,
    partnerable_agents: AgentSet,
    sex_partners: Dict,
    pwid_agents: AgentSet,
    params: ObjMap,
    rand_gen,
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
    eligible_partner_set = partnerable_agents.members - agent.partners - {agent}

    def bondtype(bond_dict):
        bonds = list(params.classes.bond_types.keys())
        probs = [bond_dict[bond].prob for bond in bonds]
        return rand_gen.choices(bonds, weights=probs, k=1).pop()

    def assort(eligible_partner_list, assort_params):
        if rand_gen.random() < assort_params.prob:
            eligible_partners = {
                partner
                for partner in eligible_partner_list
                if (
                    getattr(partner, assort_params.assort_type)
                    == assort_params.partner_type
                )
            }
        else:
            eligible_partners = {
                partner
                for partner in eligible_partner_list
                if (
                    getattr(partner, assort_params.assort_type)
                    != assort_params.partner_type
                )
            }
        return eligible_partners

    if agent.drug_use == "Inj":
        agent_bond = bondtype(params.partnership.bonds["PWID"])
    else:
        agent_bond = bondtype(params.partnership.bonds[agent.so])

    acts_allowed = params.classes.bond_types[agent_bond].acts_allowed

    if "needle" in acts_allowed:
        eligible_partner_set &= pwid_agents.members

    if "sex" in acts_allowed:
        eligible_partner_set &= sex_partners[agent.so]

    if params.features.assort_mix:
        for assort_def in params.assort_mix.values():
            if getattr(agent, assort_def.assort_type) == assort_def.agent_type:
                eligible_partner_set = assort(eligible_partner_set, assort_def)

    random_partner = utils.safe_random_choice(eligible_partner_set, rand_gen)

    return random_partner, agent_bond


@utils.memo
def sex_possible(agent_sex_type: str, partner_sex_type: str, sex_types: ObjMap) -> bool:
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
    if agent_sex_type not in sex_types:
        raise ValueError(f"Invalid agent_sex_type! {agent_sex_type}")
    if partner_sex_type not in sex_types:
        raise ValueError(f"Invalid partner_sex_type! {partner_sex_type}")

    agent_match = agent_sex_type in sex_types[partner_sex_type].sleeps_with
    partner_match = partner_sex_type in sex_types[agent_sex_type].sleeps_with

    return agent_match and partner_match


def get_partnership_duration(agent: Agent, params: ObjMap, rand_gen) -> int:
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
    dur_info = params.partnership.sex.duration

    diceroll = rand_gen.random()
    dur_bin = dur_info[5]
    for i in range(1, 5):
        if diceroll < dur_info[i].prob:
            dur_bin = dur_info[i]
            break

    duration = rand_gen.randrange(dur_bin.min, dur_bin.max, 1,)

    return duration
