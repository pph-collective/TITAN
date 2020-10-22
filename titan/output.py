#!/usr/bin/env python
# encoding: utf-8

from typing import Dict, Any, List, Iterator
from . import agent as ag
import itertools
import os

from networkx import betweenness_centrality, effective_size, density  # type: ignore
from .parse_params import ObjMap


def setup_aggregates(params: ObjMap, features, classes: List[str]) -> Dict:
    """
    Recursively create a nested dictionary of attribute values to items to count.

    Attributes are classes defined in params, the items counted are:

    * "agents"
    * "hiv"
    * "dx"
    * "aids"
    * "hiv_new"
    * "dx_new"
    * "deaths"
    * "deaths_hiv"

    Additionally, any feature enabled may have additional stats that are tracked.  See the feature's `stats` attribute.

    args:
        params: model parameters
        classes: which classes to aggregate by [params.outputs.classes]

    returns:
        dictionary of class values to counts
    """
    if classes == []:
        base_stats = {
            "agents": 0,
            "hiv": 0,
            "dx": 0,
            "aids": 0,
            "hiv_new": 0,
            "dx_new": 0,
            "deaths": 0,
            "deaths_hiv": 0,
        }

        for feature in features:
            base_stats.update({stat: 0 for stat in feature.stats})

        return base_stats

    stats = {}
    clss, *rem_clss = classes  # head, tail
    keys = [k for k in params.classes[clss]]
    for key in keys:
        stats[key] = setup_aggregates(params, features, rem_clss)

    return stats


def get_aggregates(params: ObjMap) -> Iterator:
    """
    Get iterator over all attribute combinations for output classes

    args:
        params: model parameters

    returns:
        iterator over attribute combinations
    """
    return itertools.product(
        *[list(k for k in params.classes[clss]) for clss in params.outputs.classes]
    )


def get_agg_val(stats: Dict, attrs: List, key: str) -> int:
    """
    Get the value of a key in stats given the attribute values

    args:
        stats: a nested dictionary of attributes to count
        attrs: a list of attribute values to find the count for
        key: the type of count to get the value of

    returns:
        the count of key for the given attributes
    """
    stats_item = stats
    for attr in attrs:
        stats_item = stats_item[attr]

    return stats_item[key]


def get_stats_item(stats: Dict[str, Any], attrs: List[str], agent: "ag.Agent"):
    stats_item = stats
    for attr in attrs:
        stats_item = stats_item[str(getattr(agent, attr))]

    return stats_item


def add_agent_to_stats(
    stats: Dict[str, Any], attrs: List[str], agent: "ag.Agent", key: str
):
    """
    Update the stats dictionary counts for the key given the agent's attributes

    args:
        stats: a nested dictionary of attributes to counts
        attrs: a list of attribute types (e.g. "race")
        agent: the agent whose attribute values will be evaluated
        key: the type of count to increment
    """
    stats_item = get_stats_item(stats, attrs, agent)

    stats_item[key] += 1


def get_stats(
    all_agents: "ag.AgentSet",
    deaths: List["ag.Agent"],
    params: ObjMap,
    features,
    time: int,
) -> Dict:
    """
    Get the current statistics for a model based on the population, and tracking agent sets from the model.

    args:
        all_agents: all of the agents in the population
        new_hiv_dx: agents who are newly diagnosed with hiv this timestep
        deaths: agents who died this timestep
        params: model parameters

    returns:
        nested dictionary of agent attributes to counts of various items
    """
    stats = setup_aggregates(params, features, params.outputs.classes)
    attrs = [
        clss[:-1] for clss in params.outputs.classes
    ]  # attribute version (non-plural)

    for a in all_agents:
        stats_item = get_stats_item(stats, attrs, a)

        add_agent_to_stats(stats, attrs, a, "agents")

        for feature in features:
            agent_feature = getattr(a, feature.name)
            agent_feature.set_stats(stats_item, time)

        if a.hiv:
            add_agent_to_stats(stats, attrs, a, "hiv")
            if a.hiv_time == time:
                add_agent_to_stats(stats, attrs, a, "hiv_new")
            if a.aids:
                add_agent_to_stats(stats, attrs, a, "aids")
            if a.hiv_dx:
                add_agent_to_stats(stats, attrs, a, "dx")
                if a.hiv_dx_time == time:
                    add_agent_to_stats(stats, attrs, a, "dx_new")

    for a in deaths:
        add_agent_to_stats(stats, attrs, a, "deaths")
        if a.hiv:
            add_agent_to_stats(stats, attrs, a, "deaths_hiv")

    return stats


# ================== Printer Functions =========================
# Each of the following functions takes in the time, seeds, and stats dict for that time
# and prints the appropriate stats to file


def write_report(
    file_name: str,
    run_id: str,
    t: int,
    runseed: int,
    popseed: int,
    stats: Dict,
    params: ObjMap,
    outdir: str,
):
    """
    Core function for writing reports, writes header if file is new, then data based
    on the `params` and `name_map`

    args:
        file_name: Name of the file to write, including the extension (e.g. `MyReport.txt`)
        run_id: unique identifier for this model
        t: current timestep
        runseed: integer used to seed the random number generator for the model
        popseed: integer used to seed the random number generator for the population
        stats: nested dictionary of agent attributes to counts
        params: model parameters
        outdir: path of where to save this file
    """

    def get_stat_names(stats, attrs):
        stat_ref = stats
        for i in range(len(attrs)):
            stat_ref = stat_ref[list(stat_ref.keys())[0]]

        return stat_ref

    f = open(os.path.join(outdir, file_name), "a")
    attrs = [clss[:-1] for clss in params.outputs.classes]
    stat_names = get_stat_names(stats, attrs)

    if f.tell() == 0:
        f.write("run_id\trseed\tpseed\tt\t")  # start header

        # attributes in stats
        f.write("\t".join(attrs))

        # report specific fields
        for name in stat_names:
            f.write(f"\t{name}")

        f.write("\n")

    for agg in get_aggregates(params):
        f.write(f"{run_id}\t{runseed}\t{popseed}\t{t}\t")

        f.write("\t".join(agg))  # write attribute values

        for name in stat_names:
            f.write(f"\t{(get_agg_val(stats, agg, name))}")

        f.write("\n")

    f.close()


def basicReport(
    run_id: str,
    t: int,
    runseed: int,
    popseed: int,
    stats: Dict[str, Any],
    params: ObjMap,
    outdir: str,
):
    """
    Standard report writer for basic agent statistics, columns include:

    * "agents": number of agents in the population
    * "hiv": number of agents with HIV
    * "dx": number of agents with HIV who are diagnosed
    * "aids": number of agents with AIDS
    * "hiv_new": number of agents who HIV converted this time period
    * "dx_new": number of agents with HIV who were diagnosed this time period
    * "deaths": number of agents who died this time period
    * "deaths_hiv": number of agents with HIV who died this time period

    Additionally, any feature enabled may have additional stats that are tracked.  See the feature's `stats` attribute and docs for details.
    """
    write_report("basicReport.txt", run_id, t, runseed, popseed, stats, params, outdir)


# ========================== Other Print Functions =============================


def print_components(
    run_id: str,
    t: int,
    runseed: int,
    popseed: int,
    components: List,
    outdir: str,
    races: list,
):
    """
    Write stats describing the components (sub-graphs) in a graph to file

    args:
        run_id: unique identifer for this run of the model
        t: current timestep
        runseed: integer used to seed the model's random number generator
        popseed: integer used to seed the population's random number generator
        components: a list of graph components
        outdir: path where the file should be saved
        races: the races in the population
    """
    f = open(os.path.join(outdir, f"{run_id}_componentReport_ALL.txt"), "a")

    race_count: Dict[str, int] = {}
    race_list = []
    header = ""
    for race in races:
        race_list.append(race)
        header += "\t" + race

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write(
            "run_id\trunseed\tpopseed\tt\tcompID\ttotalN\tNhiv\tNprep\tNtrtHIV"
            "\tComp_Status\tPCA\tOral\tInjectable\tAware\tnidu\tcentrality\tDensity"
            "\tEffectiveSize" + header + "\n"
        )

    comp_id = 0
    for comp in components:
        for race in races:
            race_count[race] = 0
        assert comp.number_of_nodes() >= 0
        tot_agents = (
            nhiv
        ) = ntrthiv = nprep = injectable_prep = oral = aware = pca = nidu = 0
        component_status = "control"
        trt_comp = trt_agent = False

        for agent in comp.nodes():
            tot_agents += 1
            race_count[agent.race] += 1
            if agent.hiv:
                nhiv += 1
                if agent.random_trial.treated:
                    ntrthiv += 1

            if agent.prep.active:
                nprep += 1
                if agent.prep.type == "Inj":
                    injectable_prep += 1
                elif agent.prep.type == "Oral":
                    oral += 1

            if agent.pca.suitable and agent.pca.active:
                pca += 1

            if agent.random_trial.treated:  # treatment component
                trt_agent = True

            if agent.random_trial.active:
                trt_comp = True

            if agent.pca.awareness:
                aware += 1

            if agent.drug_type == "NonInj":
                nidu += 1

        comp_centrality = (
            sum(betweenness_centrality(comp).values()) / comp.number_of_nodes()
        )
        average_size = sum(effective_size(comp).values()) / comp.number_of_nodes()
        comp_density = density(comp)

        if trt_comp:
            if trt_agent:
                component_status = "treatment"
            else:
                component_status = "treatment_no_eligible"

        race_str = ""
        for race in race_list:
            race_str += "\t" + str(race_count[race])

        f.write(
            f"{run_id}\t{runseed}\t{popseed}\t{t}\t{comp_id}\t{tot_agents}\t"
            f"{nhiv}\t"
            f"{nprep}\t{ntrthiv}\t{component_status}\t{pca}\t{oral}\t{injectable_prep}"
            f"\t{aware}\t{nidu}\t{comp_centrality:.4f}\t{comp_density:.4f}"
            f"\t{average_size:.4f}{race_str}\n"
        )

        comp_id += 1

    f.close()
