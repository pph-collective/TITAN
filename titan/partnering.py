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
        match_fns = get_match_fns(
            params.assort_mix.values(), agent, bond_type, rand_gen
        )
        # if no definitions match this agent, don't try to assort
        if len(match_fns) > 0:
            for partner in utils.safe_shuffle(eligible, rand_gen):
                if is_assortable(partner, match_fns):
                    return partner
            return None

    return utils.safe_random_choice(eligible, rand_gen)


# does an agent match the criteria of the randomly chosen assort values?
def is_assortable(agent, match_fns):
    for match_fn in match_fns:
        if not match_fn(agent):
            return False

    return True


# get the attribute as a string, recurse if nested
def get_str_attr(obj, attr):
    attrs = attr.split(".")  # returns list of all attrs
    if len(attrs) == 1:
        return str(getattr(obj, attrs[0]))
    else:
        return get_str_attr(getattr(obj, attrs.pop(0)), ".".join(attrs))


# if this assort rule applies for this agent, get the match function for a potential partner
def get_match_fns(assort_defs, agent, bond_type, rand_gen):
    match_fns = []
    for assort_def in assort_defs:
        if len(assort_def.bond_types) == 0 or bond_type in assort_def.bond_types:
            if get_partner_attr(assort_def) == "location":
                match_fns.append(get_location_match_fn(assort_def, agent, rand_gen))
            elif assort_def.agent_value == "__any__":
                match_fns.append(get_same_match_fn(assort_def, agent, rand_gen))
            elif get_str_attr(agent, assort_def.attribute) == str(
                assort_def.agent_value
            ):
                match_fns.append(get_match_fn(assort_def, rand_gen))

    return match_fns


# pick a partner type randomly given the weights
def get_partner_type(assort_def, rand_gen):
    partner_types = list(assort_def.partner_values.keys())
    partner_weights = [assort_def.partner_values[p] for p in partner_types]
    return utils.safe_random_choice(partner_types, rand_gen, weights=partner_weights)


# what partner attribute to use in assorting
def get_partner_attr(assort_def):
    if assort_def.partner_attribute == "__agent__":
        return assort_def.attribute
    else:
        return assort_def.partner_attribute


# given an assort def, randomly select the type the partner must have given the
# weights and return a function to determin if a potential partner matches it
def get_match_fn(assort_def, rand_gen):
    partner_types = list(assort_def.partner_values.keys())
    partner_type = get_partner_type(assort_def, rand_gen)
    attr = get_partner_attr(assort_def)
    if partner_type == "__other__":
        partner_types.remove("__other__")
        return lambda ag: get_str_attr(ag, attr) not in partner_types
    else:
        return lambda ag: get_str_attr(ag, attr) == partner_type


# given an assort def where agent_value == '__any__', return the match function
def get_same_match_fn(assort_def, agent, rand_gen):
    partner_type = get_partner_type(assort_def, rand_gen)
    attr = get_partner_attr(assort_def)

    agent_attribute = get_str_attr(agent, attr)
    if partner_type == "__same__":
        return lambda ag: get_str_attr(ag, attr) == agent_attribute
    elif partner_type == "__other__":
        return lambda ag: get_str_attr(ag, attr) != agent_attribute
    else:
        raise ValueError(
            "When using same-assorting, only valid partner_types are __same__ and __other__"
        )


# given an assort def where partner_attribute == 'location', return the match function
def get_location_match_fn(assort_def, agent, rand_gen):
    partner_type = get_partner_type(assort_def, rand_gen)
    attr = "location"

    agent_location = get_str_attr(agent, attr)
    agent_neighbors = agent.location.neighbors
    if assort_def.agent_value == "__any__":
        if partner_type == "__same__":
            return lambda ag: get_str_attr(ag, attr) == agent_location
        elif partner_type == "__other__":
            if "__neighbor__" in assort_def.partner_values:
                return lambda ag: get_str_attr(ag, attr) not in agent_neighbors.union(
                    [agent_location]
                )
            else:
                return lambda ag: get_str_attr(ag, attr) != agent_location
        elif partner_type == "__neighbor__":
            return lambda ag: get_str_attr(ag, attr) in agent_neighbors
        else:
            raise ValueError(
                "When using same-assorting on location, only valid partner_types are __same__, __neighbor__ and __other__"
            )
    else:
        partner_types = list(assort_def.partner_values.keys())
        if partner_type == "__other__":
            partner_types.remove("__other__")
            if "__neighbor__" in assort_def.partner_values:
                partner_types.remove("__neighbor__")
                return lambda ag: get_str_attr(ag, attr) not in agent_neighbors.union(
                    partner_types
                )
            else:
                return lambda ag: get_str_attr(ag, attr) not in partner_types
        elif partner_type == "__neighbor__":
            return lambda ag: get_str_attr(ag, attr) in agent_neighbors
        else:
            return lambda ag: get_str_attr(ag, attr) == partner_type


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
        duration = utils.safe_random_int(dur_info[i].min, dur_info[i].max, rand_gen)

    else:
        dist = params.partnership.duration[bond_type][race].distribution
        duration = int(utils.safe_dist(dist, rand_gen))

    return duration
