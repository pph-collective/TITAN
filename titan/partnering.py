#!/usr/bin/env python
# encoding: utf-8

# Imports
from typing import Optional, Tuple, Dict
from copy import copy

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

    def bondtype(bond_dict):
        bonds = list(params.classes.bond_types.keys())
        probs = [bond_dict[bond].prob for bond in bonds]
        return rand_gen.choices(bonds, weights=probs, k=1).pop()

    def assort(eligible_partners, assort_params):
        partner_types = list(assort_params.partner_values.keys())
        partner_weights = [assort_params.partner_values[p] for p in partner_types]
        partner_type = rand_gen.choices(
            partner_types, weights=partner_weights, k=1
        ).pop()

        if partner_type == "__other__":
            # remove all the specified (not selected) values and remove those
            # partners from the eligible set
            for p in partner_types:
                if p != "__other__":
                    eligible_partners = {
                        partner
                        for partner in eligible_partners
                        if str(getattr(partner, assort_params.attribute)) != p
                    }
        else:
            eligible_partners = {
                partner
                for partner in eligible_partners
                if str(getattr(partner, assort_params.attribute)) == partner_type
            }

        return eligible_partners

    eligible = copy(partnerable_agents.members)
    eligible -= agent.partners
    eligible -= {agent}

    if agent.drug_use == "Inj":
        agent_bond = bondtype(params.partnership.bonds["PWID"])
    else:
        agent_bond = bondtype(params.partnership.bonds[agent.so])

    acts_allowed = params.classes.bond_types[agent_bond].acts_allowed

    if "injection" in acts_allowed:
        eligible &= pwid_agents.members

    if "sex" in acts_allowed:
        eligible &= sex_partners[agent.so]

    if not eligible:
        return None, agent_bond

    if params.features.assort_mix:
        for assort_def in params.assort_mix.values():
            if getattr(agent, assort_def.attribute) == assort_def.agent_value:
                eligible = assort(eligible, assort_def)

    random_partner = utils.safe_random_choice(eligible, rand_gen)

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
