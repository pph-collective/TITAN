#!/usr/bin/env python
# encoding: utf-8

# Imports
from typing import Set, Optional, Tuple

from dotmap import DotMap  # type: ignore

from .agent import Agent
from .utils import safe_random_choice


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
    partner_set = set(need_partners)
    eligible_partner_set = partner_set - set(agent.partners) - {agent}
    RandomPartner: Optional[Agent]

    def bondtype(bond_dict):
        bonds = []
        probs = []
        for bond, val in bond_dict.items():
            if bond != "prob":
                bonds.append(bond)
                probs.append(val.prob)
        bonded_type_hold = rand_gen.choices(bonds, weights=probs, k=1)
        bonded_type = "".join(bonded_type_hold)

        return bonded_type

    def assort(eligible_partner_list, assort_params):
        if rand_gen.random() < assort_params["prob"]:
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

    if params.features.assort_mix:
        for assort_types in params.assort_mix.assortativity.assort_type:
            if getattr(agent, assort_types.assort_type) == assort_types["agent_type"]:
                eligible_partner_set = assort(eligible_partner_set, assort_types)

    if agent.drug_use == "Inj":
        agent_bond = bondtype(params.partnership.bonds["PWID"])
    else:
        agent_bond = bondtype(params.partnership.bonds[agent.so])

    if "needle" in params.classes.bond_types[agent_bond].acts_allowed:
        eligible_partner_set = {
            partner for partner in eligible_partner_set if partner.drug_use == "Inj"
        }

    if "sex" in params.classes.bond_types[agent_bond].acts_allowed:
        eligible_partner_set = {
            partner
            for partner in eligible_partner_set
            if sex_possible(agent.so, partner.so, params)
        }

    if "social" in params.classes.bond_types[agent_bond].acts_allowed:
        eligible_partner_set = eligible_partner_set

    if eligible_partner_set:
        random_partner = safe_random_choice(list(eligible_partner_set), rand_gen)
    else:
        random_partner = None

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
        raise ValueError("Invalid agent_sex_type! %s" % str(agent_sex_type))
    if partner_sex_type not in sex_types:
        raise ValueError("Invalid partner_sex_type! %s" % str(partner_sex_type))


    agent_match = (
        agent_sex_type in sex_types[partner_sex_type].sleeps_with
    )
    partner_match = (
        partner_sex_type in sex_types[agent_sex_type].sleeps_with
    )

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
