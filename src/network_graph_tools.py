#!/usr/bin/env python
# encoding: utf-8

import os
import random
import collections

import numpy as np
import networkx as nx
from ns.nx_agraph import graphviz_layout
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from operator import itemgetter

try:
    from .HIVABM_Population import PopulationClass
except ImportError:
    raise ImportError("Can't import PopulationClass")

from . import params


def save_adjlist(
    N_pop, graph, dir_prefix, time
):  # REVIEW there are a few places where a flag referencing this is used, the the function isn't called anywhere
    """
    :Purpose:
    Write graph G in a single-ling agdacency-list format to path
    Parameters graph:NetworkX graph
              path: string or file. Filename or file handle for data output.
              comments: marker for comment lines
              delimiter: string, optional, separator for node labels
              encoding: string, optional, text enconding
    """
    if not os.path.isdir(dir_prefix):
        os.mkdir(dir_prefix)
    OutFileName = dir_prefix + "/" + "Interactions_at_time_%d.txt" % time
    outfile = open(OutFileName, "w")
    outfile.write("Agent 1\ttime\tAgent 2\n")
    not_singletons = []

    for e in graph.edges_iter():
        not_singletons.append(e[0])
        not_singletons.append(e[1])
        outfile.write("%d\t%d\t%d\n" % (e[0], time, e[1]))
    singletons = list(set(range(N_pop)).difference(set(not_singletons)))
    for singleton in singletons:
        outfile.write("%d\t%d\n" % (singleton, time))


class NetworkClass(PopulationClass):
    def __init__(self, N, popSeed=0, netSeed=0, m_0=1, network_type="scale_free"):
        """
        :Purpose:
            This is the base class used to generate the social network
            for the other agents, i.e. . The class inherits from the PopulationClass.

        :Input:
            N : int
              Number of agents. Default: 10000

            m_0: int
              Number of nodes each node is connected to in preferential
              attachment step

            network_type: defaul is "scale_free", other options are "max_k_comp_size" and "binomial"
        """

        random.seed(netSeed)
        np.random.seed(netSeed)

        if type(N) is not int:
            raise ValueError(
                "Population size must be integer,\
                n = %s, not int"
                % (type(N))
            )
        else:
            pass
        if m_0 not in list(range(10)):
            raise ValueError("m_0 must be integer smaller than 10")
        else:
            self.m_0 = m_0
        PopulationClass.__init__(self, n=N, rSeed=popSeed)  # Create population

        self.NetworkSize = N
        if network_type == "scale_free":
            self.G = nx.Graph()
            for i in range(10):
                self.update_partner_assignments(params.PARTNERTURNOVER, self.G)
        elif network_type == "max_k_comp_size":

            def trimComponent(component, maxComponentSize):
                for ag in component.nodes:
                    if random.random() < 0.1:
                        for rel in ag._relationships:
                            # print("Removed edge:",rel)
                            rel.progress(forceKill=True)
                            self.Relationships.remove_agent(rel)
                            component.remove_edge(rel._ID1, rel._ID2)
                            self.G.remove_edge(rel._ID1, rel._ID2)

                components = sorted(
                    nx.connected_component_subgraphs(self.G), key=len, reverse=True
                )
                totNods = 0
                for comp in components:
                    cNodes = comp.number_of_nodes()
                    if cNodes > params.maxComponentSize:
                        trimComponent(comp, params.maxComponentSize)
                    elif cNodes < params.minComponentSize:
                        for a in comp.nodes():
                            try:
                                self.G.remove_node(a)
                            except:
                                pass
                    else:
                        totNods += cNodes

            self.G = nx.Graph()
            for i in range(10):
                self.update_partner_assignments(params.PARTNERTURNOVER, self.get_Graph)
            components = sorted(
                nx.connected_component_subgraphs(self.G), key=len, reverse=True
            )
            for comp in components:
                if comp.number_of_nodes() > params.maxComponentSize:
                    print("TOO BIG", comp, comp.number_of_nodes())
                    trimComponent(comp, params.maxComponentSize)
                elif comp.number_of_nodes() < params.minComponentSize:
                    print("TOO SMALL", comp, comp.number_of_nodes())
            print("Total agents in graph: ", self.G.number_of_nodes())
        else:
            print("HUIH")
            raise ValueError("Invalid network type! %s" % str(network_type))

    def write_G_edgelist(self, path):
        G = self.G
        nx.write_edgelist(G, path, data=False)

    def write_network_stats(self, t=0, path="results/network/networkStats.txt"):
        from networkx.algorithms import approximation

        G = self.G
        components = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)
        bigG = components[0]
        outfile = open(path, "w")
        outfile.write(nx.info(G))

        centDict = nx.degree_centrality(G)

        outfile.write(
            "\nNumber of connected components: {}\n".format(
                nx.number_connected_components(G)
            )
        )
        conG = nx.connected_components(G)
        tot_nodes = 0
        for t in conG:
            thisComp = len(t)
            tot_nodes += thisComp
        outfile.write(
            "Average component size: {}\n".format(
                tot_nodes * 1.0 / nx.number_connected_components(G)
            )
        )
        outfile.write("Maximum component size: {}\n".format(nx.number_of_nodes(bigG)))
        outfile.write("Degree Histogram: {}\n".format(nx.degree_histogram(G)))
        outfile.write("Graph density: {}\n".format(nx.density(G)))
        outfile.write(
            "Average node degree centrality: {}\n".format(
                sum(centDict.values()) / len(list(centDict.values()))
            )
        )

        outfile.write("Average node clustering: {}\n".format(nx.average_clustering(G)))
        outfile.close()

        comps = []
        for i in components:
            comps.append(len(i))

    def create_graph_from_agents(self, agents):
        G = self.G
        numAdded = 0
        for tmpA in agents.iter_agents():
            numAdded += 1
            G.add_node(tmpA)
        print("\tAdded %d/%d agents" % (numAdded, G.number_of_nodes()))

    def create_graph_from_relationships(self, relationships):
        G = self.G
        numAdded = 0
        for rel in relationships.iter_agents():
            numAdded += 1
            self.G.add_edge(rel._ID1, rel._ID2)
        print("Added %d/%d relationships" % (numAdded, G.number_of_edges()))

    def draw_histogram(self, t=0):
        G = self.G
        degree_sequence = sorted(
            [d for n, d in G.degree()], reverse=True
        )  # degree sequence
        degreeCount = collections.Counter(degree_sequence)
        deg, cnt = list(zip(*list(degreeCount.items())))

        fig, ax = plt.subplots()
        plt.bar(deg, cnt, width=0.80, color="b")

        plt.title("Degree Histogram\nTime: %d" % t)
        plt.ylabel("Count")
        plt.xlabel("Degree")
        plt.ylim(0, len(G.nodes))
        ax.set_xticks([d + 0.4 for d in deg])
        ax.set_xticklabels(deg)

        # draw graph in inset
        plt.axes([0.4, 0.4, 0.5, 0.5])
        Gcc = G
        pos = graphviz_layout(Gcc, prog="neato", args="")
        plt.axis("off")

        node_shape = "o"
        node_color = []
        for v in Gcc:
            if v._AIDS_bool:
                node_color.append("r")
            elif v._HIV_bool:
                node_color.append("r")
            else:
                node_color.append("g")
        nx.draw_networkx_nodes(
            Gcc, pos, node_size=20, node_color=node_color, node_shape=node_shape
        )
        nx.draw_networkx_edges(Gcc, pos, alpha=0.4)

        plt.show(block=False)
        plt.savefig("images/snapshot_%d.png" % t)
        plt.pause(0.5)
        plt.close()

    def plot_DegreeDistribution(self, time=0):
        """
        Plot the node degree distribution of the graph. \n
        INPUT: networkX graph
        """
        graph = self.G
        Gsize = graph.number_of_nodes()
        degreeList = np.array(graph.degree())
        x_degree = np.arange(max(degreeList) + 1)  # include 0 and max
        degreehist = np.bincount(degreeList)
        ix = degreehist != 0  # delete zeros
        degreehist = degreehist[ix]
        x_degree = x_degree[ix]
        plt.clf()
        plt.plot(x_degree, degreehist, "bo")
        plt.suptitle("Node degree distribution", fontsize=18)
        plt.title("Time=%d   N=%d" % (time, Gsize))
        plt.ylabel("Frequency")
        plt.xlabel("Node degree")
        plt.axis([0, 1.05 * max(degreeList), 0, 1.05 * max(degreehist)])
        plt.savefig("images/Degree_%d.png" % time)  # plt.show()

    def get_AdjacencyList(self):
        """ Return the adjacency list of the graph. """
        return self.G.adjacency_list()

    def stat_connectivity(self):
        G = self.G
        return nx.all_pairs_node_connectivity(G)

    def get_Graph(self):
        """
        Return random assortative graph produced by ``set_assortative_graph``.
        """
        return self.G

    def vizualize_network_circular(self):
        G = self.G
        pos = graphviz_layout(G, prog="twopi", args="")
        plt.figure(figsize=(8, 8))
        nx.draw(G, pos, node_size=20, alpha=0.5, node_color="blue", with_labels=False)
        plt.axis("equal")
        plt.savefig("circular_tree.png")

    def vizualize_network_spectral(self):
        G = self.G
        pos = graphviz_layout(G, prog="twopi", args="")
        plt.figure(figsize=(8, 8))
        nx.draw_spectral(G, with_labels=False)
        plt.axis("equal")
        plt.savefig("spectral_tree.png")

    def vizualize_network_graphviz(self, program="sfdp", coloring=None, time=0):
        G = self.G
        Gsize = G.number_of_nodes()

        plt.clf()
        plt.figure(figsize=(8, 8))
        plt.suptitle("Agent Network %s" % program, fontsize=18)
        plt.title("Time=%d   N=%d" % (time, Gsize))
        size = G.size()

        edge_color = "k"
        node_shape = "o"
        node_color = self.get_network_color(coloring)

        pos = graphviz_layout(G, prog="neato", args="")
        nx.draw(
            G,
            pos,
            node_size=1,
            node_color=node_color,
            node_shape=node_shape,
            edge_color=edge_color,
            with_labels=False,
            linewidths=0.25,
            width=0.25,
        )

        plt.axis("equal")
        filename = "images/Network_%d_%s_%s_%d.png" % (Gsize, program, coloring, time)
        plt.savefig(filename)

    def vizualize_network_random(self):
        G = self.G
        plt.figure(figsize=(8, 8))
        nx.draw_random(G, with_labels=False)
        plt.axis("equal")
        plt.savefig("random.png")

    def vizualize_network_ego(self):
        # Create a BA model graph
        G = self.G
        n = 1000
        m = 2
        # find node with largest degree
        node_and_degree = G.degree()
        (largest_hub, degree) = sorted(
            list(node_and_degree.items()), key=itemgetter(1)
        )[-1]
        # Create ego graph of main hub
        hub_ego = nx.ego_graph(G, largest_hub)
        # Draw graph
        nx.draw_networkx(G, with_labels=False)
        plt.savefig("ego_graph.png")

    def get_network_color(self, coloring="Sex Type"):
        G = self.G
        node_color = []
        if coloring == "SO":
            for v in G:
                tmp_sextype = v._SO
                if tmp_sextype == "HM":
                    node_color.append("b")
                elif tmp_sextype == "HF":
                    node_color.append("g")
                elif tmp_sextype == "WSW":
                    node_color.append("c")
                elif tmp_sextype == "MSM":
                    node_color.append("r")
                else:
                    raise ValueError("Check agents %s sextype %s" % (v, tmp_sextype))
        elif coloring == "DU":
            for v in G:
                tmp_drugtype = self.get_agent_characteristic(v, "Drug Type")
                if tmp_drugtype == "ND":
                    node_color.append("g")
                elif tmp_drugtype == "NIDU":
                    node_color.append("b")
                elif tmp_drugtype == "IDU":
                    node_color.append("r")
                else:
                    raise ValueError("Check agents %s drug type %s" % (v, tmp_drugtype))
        elif coloring == "Tested":
            for v in G:
                if v._HAART_bool:
                    node_color.append("g")
                elif v._tested:  # tmp_hiv == 1:
                    node_color.append("y")
                elif v._HIV_bool:  # tmp_aids == 1:
                    node_color.append("r")
                elif v._PrEP_bool:
                    node_color.append("b")
                else:
                    node_color.append("purple")
        elif coloring == "Trtmt":
            for v in G:
                if v._HIV_bool:  # tmp_aids == 1:
                    node_color.append("r")
                elif v._PrEP_bool:
                    node_color.append("g")
                elif v._treatment_bool:
                    node_color.append("y")
                else:
                    node_color.append("gray")
        elif coloring == "HIV":
            for v in G:
                if v._AIDS_bool:  # tmp_hiv == 1:
                    node_color.append("purple")
                elif v._HIV_bool:  # tmpaids == 1:
                    node_color.append("r")
                else:
                    node_color.append("g")
        elif coloring == "HR":
            for v in G:
                if v._highrisk_bool:  # tmp_hiv == 1:
                    node_color.append("r")
                elif v._everhighrisk_bool:  # tmp_aids == 1:
                    node_color.append("y")
                else:
                    node_color.append("g")
        elif coloring == "Race":
            for v in G:
                if v._race == "WHITE":
                    node_color.append("y")
                elif v._race == "BLACK":  # tmp_aids == 1:
                    node_color.append("g")
                else:
                    node_color.append("b")
        elif coloring == "MSW":
            for v in G:
                if v._race == "BLACK":
                    node_color.append("y")
                elif v._everhighrisk_bool:
                    node_color.append("b")
                elif v._race == "WHITE":
                    node_color.append("g")
                else:
                    raise ValueError("Check agents %s drug type %s" % (v, tmp_drugtype))
        else:
            raise ValueError(
                "coloring value invalid!\n%s\n \
            Only 'SO','DU','Tested', 'Trtmt', and 'HIV' allowed!"
                % str(coloring)
            )

        return node_color

    def visualize_network(
        self,
        coloring="SO",
        pos=None,
        return_layout=0,
        node_size=None,
        iterations=1,
        curtime=0,
        txtboxLabel=0,
        label="Network",
    ):
        """
        :Purpose:
            Visualize the network using the spring layout (default). \n

        :Input:
            graph : networkX graph
        """
        G = self.G
        print(("\tPlotting {} colored by {}...").format(label, coloring))
        fig = plt.figure()
        ax = fig.add_axes([0, 0, 1, 1])
        fig.clf()

        # build a rectangle in axes coords
        left, width = 0.0, 1.0
        bottom, height = 0.0, 1.0
        right = left + width
        top = bottom + height

        fig = plt.figure()
        ax = fig.add_axes([0, 0, 1, 1])

        # axes coordinates are 0,0 is bottom left and 1,1 is upper right
        p = patches.Rectangle(
            (left, bottom),
            width,
            height,
            fill=False,
            transform=ax.transAxes,
            clip_on=False,
        )

        ax.add_patch(p)

        if not pos:
            pos = graphviz_layout(G, prog="neato", args="")

        edge_color = "k"
        node_shape = "o"

        # node color to by type
        node_color = self.get_network_color(coloring)

        # node size indicating node degree
        NodeSize = []
        if node_size:
            for v in G:
                NodeSize.append(node_size)
        else:
            for v in G:
                NodeSize.append((10 * G.degree(v)) ** (1.0))

        # draw:
        drawMe = nx.draw(
            G,
            pos,
            node_size=NodeSize,
            node_color=node_color,
            node_shape=node_shape,
            edge_color=edge_color,
            with_labels=False,
            linewidths=0.5,
            width=0.5,
        )

        textstr = "\n".join(
            (r"N infection=%.2f" % (txtboxLabel,), r"Time=%.2f" % (curtime,))
        )

        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle="round", facecolor="wheat", alpha=0.9)

        # place a text box in upper right in axes coords
        ax.text(
            right - 0.025,
            top - 0.025,
            textstr,
            horizontalalignment="right",
            verticalalignment="top",
            transform=ax.transAxes,
            bbox=props,
        )

        filename = "images/%s_%d_%s_%d.png" % (
            label,
            G.number_of_nodes(),
            coloring,
            curtime,
        )

        fig.savefig(filename)

        if return_layout:
            return pos
