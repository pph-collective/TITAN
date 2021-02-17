#!/usr/bin/env python
# encoding: utf-8

# Imports
from typing import Optional, Dict, Set
from copy import copy

import numpy as np  # type: ignore

from . import agent  # import Agent, AgentSet
from . import utils
from . import parse_params


def select_partner(
    agent: "agent.Agent",
    partnerable_agents: Set["agent.Agent"],
    sex_partners: Dict,
    pwid_agents: "agent.AgentSet",
    params: "parse_params.ObjMap",
    rand_gen,
    bond_type: str,
) -> Optional["agent.Agent"]:
    """
    Get a partner for the agent.

    args:
        agent : agent in need of a partner
        partnerable_agents: agents that can be selected as a partner
        sex_partners: mapping from sex_type to agents in the population that can sleep with that sex_type
        pwid_agents: agents with `drug_type==="Inj"`
        params: model parameters
        rand_gen: random number generator
        bond_type: type of relationship that is being formed with the partner

    returns:
        new partner or `None`
    """
    eligible = copy(partnerable_agents)
    eligible -= agent.get_partners()
    eligible -= {agent}

    acts_allowed = params.classes.bond_types[bond_type].acts_allowed

    if "injection" in acts_allowed:
        eligible &= pwid_agents.members

    if "sex" in acts_allowed:
        eligible &= sex_partners[agent.sex_type]

    # short circuit to avoid attempting to assort with no eligible partners
    if not eligible:
        return None

    if params.features.assort_mix:
        random_partner = None
        assort_attrs = get_assort_attrs(params.assort_mix.values(), agent, rand_gen)
        for partner in utils.safe_shuffle(eligible, rand_gen):
            if is_assortable(partner, assort_attrs):
                random_partner = partner
                break
    else:
        random_partner = utils.safe_random_choice(eligible, rand_gen)

    return random_partner


def is_assortable(agent, assort_attrs):
    for attr, defn in assort_attrs.items():
        if not defn["compare"](agent, attr, defn["value"]):
            return False

    return True


def get_assort_attrs(assort_defs, agent, rand_gen):
    assort_attrs = {}
    for assort_def in assort_defs:
        if getattr(agent, assort_def.attribute) == assort_def.agent_value:
            assort_attrs[assort_def.attribute] = get_assort_attr_value(
                assort_def, rand_gen
            )

    return assort_attrs


def get_assort_attr_value(assort_def, rand_gen):
    partner_types = list(assort_def.partner_values.keys())
    partner_weights = [assort_def.partner_values[p] for p in partner_types]
    partner_type = utils.safe_random_choice(
        partner_types, rand_gen, weights=partner_weights
    )
    if partner_type == "__other__":
        compare = lambda ag, attr, v: str(getattr(ag, attr)) not in v
        partner_types.remove("__other__")
        val = partner_types
    else:
        compare = lambda ag, attr, v: str(getattr(ag, attr)) == v
        val = partner_type

    return {"compare": compare, "value": val}


@utils.memo
def sex_possible(
    agent_sex_type: str, partner_sex_type: str, sex_types: "parse_params.ObjMap"
) -> bool:
    """
    Determine if sex is possible.

    args:
        agent_sex_type: name of agent's sex type
        partner_sex_type: name of partner's sex type
        sex_types: params defining in scope sex types

    returns:
        whether sex is possible between agents of these sex types
    """

    # Check input
    if agent_sex_type not in sex_types:
        raise ValueError(f"Invalid agent_sex_type! {agent_sex_type}")
    if partner_sex_type not in sex_types:
        raise ValueError(f"Invalid partner_sex_type! {partner_sex_type}")

    agent_match = agent_sex_type in sex_types[partner_sex_type].sleeps_with
    partner_match = partner_sex_type in sex_types[agent_sex_type].sleeps_with

    return agent_match and partner_match


def get_mean_rel_duration(params: "parse_params.ObjMap"):
    """
    Find the average partnership duration by bond type

    args:
        params: The current model's parameters
    """
    mean_rel_duration: Dict[str, Dict] = {}
    for bond in params.partnership.duration:
        mean_rel_duration[bond] = {}
        for race in params.classes.races:
            if params.partnership.duration[bond][race].type == "bins":
                weights = []
                vals = []
                dur_bins = params.partnership.duration[bond][race].bins
                for bins in dur_bins:
                    if bins > 1:
                        weights.append(dur_bins[bins].prob - dur_bins[bins - 1].prob)
                    else:
                        weights.append(dur_bins[bins].prob)
                    vals.append(np.average([dur_bins[bins].min, dur_bins[bins].max]))
                mean_rel_duration[bond][race] = np.average(vals, weights=weights)
            else:
                mean_rel_duration[bond][race] = params.partnership.duration[bond][
                    race
                ].distribution.mean
            assert (
                mean_rel_duration[bond][race] > 0
            ), "All bonds must have a positive duration!"

    return mean_rel_duration


def get_partnership_duration(
    params: "parse_params.ObjMap", rand_gen, bond_type: str, race: Optional[str]
) -> int:
    """
    Get duration of a relationship drawn from bins or a distribution per the params [params.partnership.duration]

    args:
        params: model parameters
        rand_gen: np random number generator
        bond_type: type of bond for the relationship whose duration is being determined

    returns:
        number of time steps the partnership should endure
    """

    if params.partnership.duration[bond_type][race].type == "bins":
        dur_info = params.partnership.duration[bond_type][race].bins
        i = utils.get_independent_bin(rand_gen, dur_info)
        duration = utils.safe_rand_int(dur_info[i].min, dur_info[i].max, rand_gen)

    else:
        dist = params.partnership.duration[bond_type][race].distribution
        duration = int(utils.safe_dist(dist, rand_gen))

    return duration
