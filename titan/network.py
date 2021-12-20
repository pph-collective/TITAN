#!/usr/bin/env python
# encoding: utf-8

import os
from typing import List, Optional

import networkx as nx  # type: ignore
from networkx.drawing.nx_agraph import graphviz_layout  # type: ignore


class NetworkGraphUtils:
    def __init__(self, graph: nx.Graph):
        """
        This is the base class used to gather statistics from an
            exsting networkx graph object.

        args:
            graph : NetworkX graph object (typically attached to a Population's `graph` attribute)
        """
        self.G = graph

    def connected_components(self) -> List:
        """
        Get connected components in the graph

        returns:
            list of connected components
        """
        return list(self.G.subgraph(c).copy() for c in nx.connected_components(self.G))

    def write_graph_edgelist(self, path: str, id, time):
        """
        Writes a pipe-delimited edge list to the file `<id>_Edgelist_t<time>.txt`

        args:
            path: directory where the file should be saved
            id: identifier for the network, typically the model's `id`
            time: timestep the edgelist is being written at
        """
        file_path = os.path.join(path, f"{id}_Edgelist_t{time}.txt")
        # Write edgelist with bond type
        nx.write_edgelist(self.G, file_path, delimiter="|", data=["type"])

    def write_network_stats(self, path: str, id, time):
        """
        Writes network statistics to the file `<id>_NetworkStats_t<time>.txt`

        args:
            path: directory where the file should be saved
            id: identifier for the network, typically the model's `id`
            time: timestep the edgelist is being written at
        """
        file_path = os.path.join(path, f"{id}_NetworkStats_t{time}.txt")

        components = sorted(self.connected_components(), key=len, reverse=True)

        outfile = open(file_path, "w")
        outfile.write(nx.info(self.G))

        cent_dict = nx.degree_centrality(self.G)

        outfile.write(
            "\nNumber of connected components: {}\n".format(
                nx.number_connected_components(self.G)
            )
        )

        tot_nodes = 0
        for c in components:
            tot_nodes += c.number_of_nodes()

        outfile.write(
            "Average component size: {}\n".format(
                tot_nodes * 1.0 / nx.number_connected_components(self.G)
            )
        )
        outfile.write(
            "Maximum component size: {}\n".format(nx.number_of_nodes(components[0]))
        )
        outfile.write("Degree Histogram: {}\n".format(nx.degree_histogram(self.G)))
        outfile.write("Graph density: {}\n".format(nx.density(self.G)))
        outfile.write(
            "Average node degree centrality: {}\n".format(
                sum(cent_dict.values()) / len(list(cent_dict.values()))
            )
        )

        outfile.write(
            "Average node clustering: {}\n".format(nx.average_clustering(self.G))
        )
        outfile.close()

    def get_network_color(self, coloring) -> List[str]:
        """
        Get a vector of node colors for plotting this graph based on the type of coloring requested

        args:
            coloring: attribute to color nodes by (e.g. "race", "diagnosed")

        returns:
            list of colors in node order
        """
        G = self.G
        node_color = []

        # attribute based coloring
        color_order = ["b", "g", "c", "r", "y", "purple", "gray"]
        if hasattr(list(G.nodes)[0], coloring):
            attrs: List[str] = []
            for v in G:
                val = getattr(v, coloring)
                if val not in attrs:
                    attrs.append(val)
                node_color.append(color_order[attrs.index(val)])

            return node_color

        # hard coded coloring schemes
        if coloring == "Diagnosed":
            for v in G:
                if v.haart:
                    node_color.append("g")
                elif v.hiv_dx:  # tmp_hiv == 1:
                    node_color.append("y")
                elif v.hiv:  # tmp_aids == 1:
                    node_color.append("r")
                elif v.prep:
                    node_color.append("b")
                else:
                    node_color.append("purple")
        elif coloring == "Trtmt":
            for v in G:
                if v.hiv:  # tmp_aids == 1:
                    node_color.append("r")
                elif v.prep:
                    node_color.append("g")
                elif v.intervention_ever:
                    node_color.append("y")
                else:
                    node_color.append("gray")
        elif coloring == "HIV":
            for v in G:
                if v.aids:  # tmp_hiv == 1:
                    node_color.append("purple")
                elif v.hiv:  # tmpaids == 1:
                    node_color.append("r")
                else:
                    node_color.append("g")
        elif coloring == "HR":
            for v in G:
                if v.high_risk:  # tmp_hiv == 1:
                    node_color.append("r")
                elif v.high_risk_ever:  # tmp_aids == 1:
                    node_color.append("y")
                else:
                    node_color.append("g")
        else:
            raise ValueError(
                "coloring value invalid!\n{coloring}\n \
            Only 'Diagnosed', 'Trtmt', 'HR', 'HIV', or an Agent attribute allowed!"
            )

        return node_color

