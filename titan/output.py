#!/usr/bin/env python
# encoding: utf-8

from typing import Dict, Any, List, Iterator
import itertools
import os

import networkx as nx  # type: ignore
from numpy import mean  # type: ignore

from .parse_params import ObjMap
from . import utils
from . import agent as ag


def setup_aggregates(
    params: ObjMap, reportables, classes: List[str], exits: Dict[str, List["ag.Agent"]]
) -> Dict:
    """
    Recursively create a nested dictionary of attribute values to items to count.

    Attributes are classes defined in params, the items counted are:

    * "agents"
    * "[exit_class]_hiv"
    * exit by exit class

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
        }
        for exit in params.classes.exit.keys():
            base_stats[exit] = 0
            base_stats[exit + "_hiv"] = 0

        for reportable in reportables:
            base_stats.update({stat: 0 for stat in reportable.stats})

        return base_stats

    stats = {}
    clss, *rem_clss = classes  # head, tail
    keys = [k for k in params.classes[clss]]

    for key in keys:
        stats[key] = setup_aggregates(params, reportables, rem_clss, exits)

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
    """
    Get the leaf node of the stats dictionary for the given attributes and agent.

    args:
        stats: a nested dictionary of attributes to count
        attrs: a list of attribute values to find the count for
        agent: The agent to get the leaf node for

    returns:
        a stats_item dictionary of keys to counts
    """
    stats_item = stats
    for attr in attrs:
        stats_item = stats_item[str(getattr(agent, attr))]

    return stats_item


def add_agent_to_stats(stats_item: Dict[str, int], key: str):
    """
    Update the stats dictionary counts for the key given the agent's attributes

    args:
        stats_item: the leaf node of a nested dictionary of attributes to counts
        key: the type of count to increment
    """
    stats_item[key] += 1


def get_stats(
    all_agents: "ag.AgentSet",
    exits: Dict[str, List["ag.Agent"]],
    params: ObjMap,
    exposures,
    features,
    time: int,
) -> Dict:
    """
    Get the current statistics for a model based on the population, and tracking agent sets from the model.

    args:
        all_agents: all of the agents in the population
        new_hiv.dx: agents who are newly diagnosed with hiv this timestep
        exits: dictionary indicated how many agents exited by exit class this timestep
        params: model parameters

    returns:
        nested dictionary of agent attributes to counts of various items
    """
    reportables = exposures + features
    stats = setup_aggregates(params, reportables, params.outputs.classes, exits)

    # attribute names (non-plural)
    attrs = [clss[:-1] for clss in params.outputs.classes]

    for a in all_agents:
        stats_item = get_stats_item(stats, attrs, a)

        add_agent_to_stats(stats_item, "agents")

        for reportable in reportables:
            agent_feature = getattr(a, reportable.name)
            agent_feature.set_stats(stats_item, time)

    for exit in exits:
        for a in exits[exit]:
            stats_item = get_stats_item(stats, attrs, a)
            add_agent_to_stats(stats_item, exit)
            if a.hiv.active:  # type: ignore[attr-defined]
                add_agent_to_stats(stats_item, exit + "_hiv")

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
        # don't write row if no agents are in it
        if get_agg_val(stats, agg, "agents") > 0:
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
    * "[exit]": number of agents who exited this time period by exit class
    * "[exit]_hiv": number of agents with HIV who exited this time period by exit class

    Additionally, any feature enabled may have additional stats that are tracked.  See the feature's `stats` attribute and docs for details.
    """
    write_report("basicReport.txt", run_id, t, runseed, popseed, stats, params, outdir)


# ========================== Other Print Functions =============================


def effective_size(comp):
    total_size = 0
    for node in comp.nodes():
        num_neighbors = len(comp[node])
        if num_neighbors == 0:
            continue

        E = comp.subgraph(comp.neighbors(node))
        total_size += num_neighbors - (2 * E.size()) / num_neighbors

    return total_size


def print_components(
    run_id: str,
    t: int,
    runseed: int,
    popseed: int,
    components: List,
    outdir: str,
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

    # if this is a new file, write the header info
    if f.tell() == 0:
        f.write(
            "run_id\trunseed\tpopseed\tt\tcomponent"
            "\tdensity\tEffectiveSize\tdeg_cent\n"
        )

    for (id, comp) in enumerate(components):
        num_nodes = comp.number_of_nodes()

        average_size = effective_size(comp) / num_nodes
        comp_density = nx.density(comp)

        deg_cent = mean(list(nx.degree_centrality(comp).values()))

        f.write(
            f"{run_id}\t{runseed}\t{popseed}\t{t}\t{id}\t"
            f"\t{comp_density:.4f}"
            f"\t{average_size:.4f}\t{deg_cent}\n"
        )

    f.close()


def write_graph_edgelist(graph, path: str, id, time):
    """
    Writes a pipe-delimited edge list to the file `<id>_Edgelist_t<time>.txt`

    args:
        path: directory where the file should be saved
        id: identifier for the network, typically the model's `id`
        time: timestep the edgelist is being written at
    """
    file_path = os.path.join(path, f"{id}_Edgelist_t{time}.txt")
    # Write edgelist with bond type
    nx.write_edgelist(graph, file_path, delimiter="|", data=["type"])


def write_network_stats(graph, path: str, id, time):
    """
    Writes network statistics to the file `<id>_NetworkStats_t<time>.txt`

    args:
        path: directory where the file should be saved
        id: identifier for the network, typically the model's `id`
        time: timestep the edgelist is being written at
    """
    file_path = os.path.join(path, f"{id}_NetworkStats_t{time}.txt")

    components = utils.connected_components(graph)

    outfile = open(file_path, "w")
    outfile.write(nx.info(graph))

    cent_dict = nx.degree_centrality(graph)

    outfile.write(
        "\nNumber of connected components: {}\n".format(
            nx.number_connected_components(graph)
        )
    )

    tot_nodes = 0
    for c in components:
        tot_nodes += c.number_of_nodes()

    outfile.write(f"Number of nodes: {tot_nodes}\n")

    outfile.write(
        "Average component size: {}\n".format(
            tot_nodes * 1.0 / nx.number_connected_components(graph)
        )
    )
    outfile.write(
        "Maximum component size: {}\n".format(nx.number_of_nodes(components[0]))
    )
    outfile.write("Degree Histogram: {}\n".format(nx.degree_histogram(graph)))
    outfile.write("Graph density: {}\n".format(nx.density(graph)))
    outfile.write(
        "Average node degree centrality: {}\n".format(
            sum(cent_dict.values()) / len(list(cent_dict.values()))
        )
    )

    outfile.write("Average node clustering: {}\n".format(nx.average_clustering(graph)))
    outfile.close()
