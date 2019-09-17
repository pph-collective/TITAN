#!/usr/bin/env python
# encoding: utf-8

"""
*****************************************************************************
Author(s):	Maximilian King  (previous authors: Lars Seemann - lseemann@uh.edu)
Email: Maximilian_King@brown.edu
Organization: Marshall Lab, Department of Epidemiology - Brown University

Description:
    Module responsible for ABM simulation events. Operates main loop over simulation run.
    Handles agent pairing, interaction, disease propagation, interventions, deaths, etc.


Copyright (c) 2016, Maximilian King
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************
"""
__author__ = "Lars Seemann (lseemann@uh.edu)"

import os
import random
import collections
import itertools
import unittest

# import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
from operator import itemgetter

# try:
#     import pygraphviz
#     from networkx.drawing.nx_agraph import graphviz_layout
# except ImportError:
#     try:
#         import pydotplus
#         from networkx.drawing.nx_pydot import graphviz_layout
#     except ImportError:
#         raise ImportError("This example needs Graphviz and either "
#                           "PyGraphviz or PyDotPlus")

try:
    from .HIVABM_Population import PopulationClass
except ImportError:
    raise ImportError("Can't import PopulationClass")

try:
    from .agent import *
except ImportError:
    raise ImportError("Can't import agent")

try:
    from .ABM_partnering import *
except ImportError as e:
    raise ImportError("Can't import ABM_partnering! %s" % str(e))

# def save_adjlist(graph, dir_prefix, time):
def save_adjlist(N_pop, graph, dir_prefix, time):
    """
    :Purpose:
    Write graph G in a single-ling agdacency-list format to path
    Parameters graph:NetworkX graph
              path: string or file. Filename or file handle for data output.
              comments: marker for comment lines
              delimiter: string, optional, separator for node labels
              encoding: string, optional, text enconding
    """
    # nx.write_adjlist(graph,(dir_prefix +'/'+'test_adjlist_at_time_%d.txt'%time))
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


def _random_subset(seq, m):
    """ Return m unique elements from seq.

    This differs from random.sample which can return repeated
    elements if seq holds repeated elements.
    """
    targets = set()
    while len(targets) < m:
        x = random.choice(seq)
        targets.add(x)
    return targets


def my_erdos_renyi_binomial_random_graph(node_list=None, p=0.0, seed=None, directed=False):

    """Return a random graph G_{n,p} (Erdős-Rényi graph, binomial graph).

    Chooses each of the possible edges with probability p.

    This is also called binomial_graph and erdos_renyi_graph.

    Parameters
    ----------
    node_list : list of nodes, optional
        The nodes used to build the graph.
    p : float
        Probability for edge creation.
    seed : int, optional
        Seed for random number generator (default=None).
    directed : bool, optional (default=False)
        If True return a directed graph

    See Also
    --------
    fast_gnp_random_graph

    Notes
    -----
    This is an O(n^2) algorithm.  For sparse graphs (small p) see
    fast_gnp_random_graph for a faster algorithm.

    References
    ----------
    .. [1] P. Erdős and A. Rényi, On Random Graphs, Publ. Math. 6, 290 (1959).
    .. [2] E. N. Gilbert, Random Graphs, Ann. Math. Stat., 30, 1141 (1959).
    """
    if directed:
        G = nx.DiGraph()
    else:
        G = nx.Graph()
    n = len(node_list)
    G.add_nodes_from(node_list)
    G.name = "my_erdos_renyi_binomial_random_graph(%s,%s)" % (n, p)
    if p <= 0:
        return G
    if p >= 1:
        return complete_graph(n, create_using=G)

    if not seed is None:
        random.seed(seed)

    if G.is_directed():
        edges = itertools.permutations(node_list, 2)
    else:
        edges = itertools.combinations(node_list, 2)

    for e in edges:
        if random.random() < p:
            G.add_edge(*e)
    return G


def my_barabasi_albert_graph(n, m, node_list=None, seed=None):
    """Return random graph using Barabási-Albert preferential attachment model.

    Modified by Lars Seemann (lseemann@uh.edu).

    A graph of n nodes is grown by attaching new nodes each with m
    edges that are preferentially attached to existing nodes with high
    degree.

    Parameters
    ----------
    n : int
        Number of nodes
    m : int
        Number of edges to attach from a new node to existing nodes
    node_list : list of nodes, optional
        The nodes used to build the graph.
    seed : int, optional
        Seed for random number generator (default=None).

    Returns
    -------
    G : Graph

    Notes
    -----
    The initialization is a graph with with m nodes and no edges.

    References
    ----------
    .. [1] A. L. Barabási and R. Albert "Emergence of scaling in
       random networks", Science 286, pp 509-512, 1999.
    """

    if m < 1 or m >= n:
        raise nx.NetworkXError("Barabási-Albert network must have m>=1 and m<n, m=%d,n=%d" % (m, n))

    if node_list is not None:
        if n != len(node_list):
            raise nx.NetworkXError(
                "node_list must have smae length as n! len(node_list)=%d,n=%d" % (len(node_list), n)
            )
        else:
            random.shuffle(node_list)  # Shuffle node_list
            targets = node_list[:m]
    else:
        targets = random.shuffle(list(range(m)))  # Target nodes for new edges

    if seed is not None:
        random.seed(seed)

    G = nx.Graph()  # networkX graph
    G.add_nodes_from(tavideorgets[:m])  # Add m initial nodes (m0 in barabasi-speak)
    G.name = "mod_barabasi_albert_graph(%s,%s)" % (n, m)
    repeated_nodes = []  # List of existing nodes, with nodes repeated once for each adjacent edge
    Num_source = m  # Start adding the other n-m nodes. The first node is m.
    while Num_source < n:
        source = node_list[Num_source]
        G.add_edges_from(list(zip([source] * m, targets)))  # Add edges to m nodes from the source.
        repeated_nodes.extend(targets)  # Add one node to the list for each new edge just created.
        repeated_nodes.extend(
            [source] * m
        )  # And the new node "source" has m edges to add to the list.
        # Now choose m unique nodes from the existing nodes
        # Pick uniformly from repeated_nodes (preferential attachement)
        targets = _random_subset(repeated_nodes, m)
        Num_source += 1
    return G


class populationGraph:
    def create_graph_from_relationships(self, relationships):
        G = nx.Graph()

        for rel in relationships.iter_agents():
            G.add_edge(rel._ID1, rel._ID2)

        return G


class NetworkClass(PopulationClass):
    def __init__(self, N, popSeed=0, netSeed=0, m_0=1, network_type="scale_free", node_list=None):
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
                update_partner_assignments(self, params.PARTNERTURNOVER, self.get_Graph)
            # scale free Albert Barabsai Graph
            # self.G = my_barabasi_albert_graph(self.NetworkSize,m_0, node_list = node_list)
        elif network_type == "max_k_comp_size":

            def trimComponent(component, maxComponentSize):
                # print "TRIMMING", component.number_of_nodes(), component.number_of_edges()
                for ag in component.nodes:
                    if random.random() < 0.1:
                        for rel in ag._relationships:
                            # print("Removed edge:",rel)
                            rel.progress(forceKill=True)
                            self.Relationships.remove_agent(rel)
                            component.remove_edge(rel._ID1, rel._ID2)
                            self.G.remove_edge(rel._ID1, rel._ID2)
                            # print(self.G.number_of_nodes())

                # print "RESULT:", component.number_of_nodes(), component.number_of_edges()

                components = sorted(nx.connected_component_subgraphs(self.G), key=len, reverse=True)
                totNods = 0
                for comp in components:
                    cNodes = comp.number_of_nodes()
                    if cNodes > params.maxComponentSize:
                        # print "TOO BIG", comp, comp.number_of_nodes()
                        trimComponent(comp, params.maxComponentSize)
                    elif cNodes < params.minComponentSize:
                        # print "TOO SMALL", comp.nodes(), comp.number_of_nodes()
                        for a in comp.nodes():
                            try:
                                self.G.remove_node(a)
                            except:
                                pass
                    else:
                        totNods += cNodes
                # print "Total agents in graph: ",totNods,self.G.number_of_nodes()

            self.G = nx.Graph()
            for i in range(10):
                update_partner_assignments(self, params.PARTNERTURNOVER, self.get_Graph)
            components = sorted(nx.connected_component_subgraphs(self.G), key=len, reverse=True)
            for comp in components:
                if comp.number_of_nodes() > params.maxComponentSize:
                    print("TOO BIG", comp, comp.number_of_nodes())
                    trimComponent(comp, params.maxComponentSize)
                elif comp.number_of_nodes() < params.minComponentSize:
                    print("TOO SMALL", comp, comp.number_of_nodes())
            print("Total agents in graph: ", self.G.number_of_nodes())
        elif network_type == "binomial":
            self.G = my_erdos_renyi_binomial_random_graph(
                node_list=self.NormalAgents, p=0.5, seed=None, directed=False
            )
        else:
            print("HUIH")
            raise ValueError("Invalid network type! %s" % str(network_type))

    def write_G_edgelist(self, path):
        G = self.G
        nx.write_edgelist(G, path, data=False)

    def write_network_stats(self, t=0):
        from networkx.algorithms import approximation

        G = self.G
        components = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)
        bigG = components[0]
        outfile = open("results/network/networkStats.txt", "w")
        outfile.write(nx.info(G))

        centDict = nx.degree_centrality(G)

        outfile.write(
            "\nNumber of connected components: {}\n".format(nx.number_connected_components(G))
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
        # outfile.write("Average node connectivity: {}\n".format(nx.node_connectivity(G)))
        outfile.write("Average node clustering: {}\n".format(nx.average_clustering(G)))
        outfile.close()

        # # Betweenness centrality
        # bet_cen = nx.betweenness_centrality(bigG)
        # # Closeness centrality
        # clo_cen = nx.closeness_centrality(bigG)
        # # Eigenvector centrality
        # #eig_cen = nx.eigenvector_centrality(bigG)
        # print bet_cen
        # print clo_cen

        comps = []
        for i in components:
            comps.append(len(i))

        # print np.histogram(comps)

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
        degree_sequence = sorted([d for n, d in G.degree()], reverse=True)  # degree sequence
        # print "Degree sequence", degree_sequence
        degreeCount = collections.Counter(degree_sequence)
        deg, cnt = list(zip(*list(degreeCount.items())))

        fig, ax = plt.subplots()
        plt.bar(deg, cnt, width=0.80, color="b")

        # degree_sequence = sorted([d for n, d in G.degree() if n._HIV_bool], reverse=True)  # degree sequence
        # # print "Degree sequence", degree_sequence
        # degreeCount = collections.Counter(degree_sequence)
        # deg, cnt = zip(*degreeCount.items())
        # print deg
        # plt.bar(deg,cnt,color="r")
        # cntHIV =

        plt.title("Degree Histogram\nTime: %d" % t)
        plt.ylabel("Count")
        plt.xlabel("Degree")
        plt.ylim(0, len(G.nodes))
        ax.set_xticks([d + 0.4 for d in deg])
        ax.set_xticklabels(deg)

        # draw graph in inset
        plt.axes([0.4, 0.4, 0.5, 0.5])
        Gcc = G  # sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)[0]
        # pos = nx.spring_layout(G, iterations=10)
        # pos = nx.circular_layout(G)
        # pos = nx.shell_layout(G)
        # ‘neato’|’dot’|’twopi’|’circo’|’fdp’|’nop’
        pos = graphviz_layout(Gcc, prog="neato", args="")
        plt.axis("off")

        node_shape = "o"
        node_color = []
        for v in Gcc:
            # tmp_aids = self.get_agent_characteristic(v, 'AIDS')
            # tmp_hiv = self.get_agent_characteristic(v, 'HIV')
            if v._AIDS_bool:  # tmp_hiv == 1:
                node_color.append("r")
            elif v._HIV_bool:  # tmp_aids == 1:
                node_color.append("r")
            else:
                node_color.append("g")
        # nx.draw(G, node_size=5, node_color=node_color)
        nx.draw_networkx_nodes(Gcc, pos, node_size=20, node_color=node_color, node_shape=node_shape)
        nx.draw_networkx_edges(Gcc, pos, alpha=0.4)
        # nx.draw_networkx_labels(G, pos, t, font_size=16)

        # line_ani = animation.FuncAnimation(fig, interval=50, blit=True)
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
        # print("\tNetwork size = "+str(self.NetworkSize))
        degreeList = np.array(graph.degree())
        x_degree = np.arange(max(degreeList) + 1)  # include 0 and max
        degreehist = np.bincount(degreeList)
        ix = degreehist != 0  # delete zeros
        degreehist = degreehist[ix]
        x_degree = x_degree[ix]
        # plt.ion()
        plt.clf()
        # print node degree distribution
        # plt.loglog(x_degree, degreehist, 'bo')
        # plt.semilogy(x_degree, degreehist, 'bo')
        plt.plot(x_degree, degreehist, "bo")
        plt.suptitle("Node degree distribution", fontsize=18)
        plt.title("Time=%d   N=%d" % (time, Gsize))
        plt.ylabel("Frequency")
        plt.xlabel("Node degree")
        # plt.axis('equal')
        # plt.axis([min(degreeList)-(0.05*max(degreeList)),
        #          max(degreeList)+(0.05*max(degreeList)),
        #          0.9, max(degreehist)+(0.05* max(degreehist))])
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
        # plt.show()

    def vizualize_network_spectral(self):
        G = self.G
        pos = graphviz_layout(G, prog="twopi", args="")
        plt.figure(figsize=(8, 8))
        nx.draw_spectral(G, with_labels=False)
        plt.axis("equal")
        plt.savefig("spectral_tree.png")
        # plt.show()

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
        # node_color = []
        # if coloring == 'Tested':
        #     for v in G:
        #         if v._tested:#tmp_hiv == 1:
        #             node_color.append('r')
        #         elif v._HIV_bool: #tmp_aids == 1:
        #             node_color.append('r')
        #         else:
        #             node_color.append('g')

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

        # nx.draw_graphviz(G,
        #                  node_color=node_color,
        #                  node_size=int(10000.0/size),
        #                  alpha=1.0,
        #                  prog=program,
        #                  with_labels=False,
        #                  linewidths = 0.25,
        #                  width = 0.25)
        plt.axis("equal")
        filename = "images/Network_%d_%s_%s_%d.png" % (Gsize, program, coloring, time)
        plt.savefig(filename)
        # plt.show()

    def vizualize_network_random(self):
        G = self.G
        plt.figure(figsize=(8, 8))
        nx.draw_random(G, with_labels=False)
        plt.axis("equal")
        plt.savefig("random.png")
        # plt.show()

    def vizualize_network_ego(self):
        # Create a BA model graph
        G = self.G
        n = 1000
        m = 2
        # G = nx.generators.barabasi_albert_graph(n, m)
        # find node with largest degree
        node_and_degree = G.degree()
        (largest_hub, degree) = sorted(list(node_and_degree.items()), key=itemgetter(1))[-1]
        # Create ego graph of main hub
        hub_ego = nx.ego_graph(G, largest_hub)
        # Draw graph
        nx.draw_networkx(G, with_labels=False)
        # pos = nx.spring_layout(hub_ego)
        # nx.draw(hub_ego, pos, node_color='b', node_size=50, with_labels=False)
        # Draw ego as large and red
        # nx.draw_networkx_nodes(hub_ego, pos, nodelist=[largest_hub], node_size=300, node_color='r')
        plt.savefig("ego_graph.png")
        # plt.show()

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
                # tmp_aids = self.get_agent_characteristic(v, 'AIDS')
                # tmp_hiv = self.get_agent_characteristic(v, 'HIV')
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
            (left, bottom), width, height, fill=False, transform=ax.transAxes, clip_on=False
        )

        ax.add_patch(p)
        # plt.figure(figsize=(8,8))

        if not pos:
            # pos=nx.spring_layout(G,iterations=iterations)
            # pos = nx.circular_layout(G)
            # pos=nx.shell_layout(G)
            # pos=nx.random_layout(G)
            # pos=nx.spectral_layout(G)
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

        # nx.draw_networkx_nodes(self.graph,pos,node_size=NodeSize)
        # nx.draw_networkx_edges(self.graph,pos,alpha=0.4)
        # print G.number_of_nodes()
        # print curtime

        textstr = "\n".join((r"N infection=%.2f" % (txtboxLabel,), r"Time=%.2f" % (curtime,)))

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

        # plt.axis('equal')
        # plt.axis('off')
        filename = "images/%s_%d_%s_%d.png" % (label, G.number_of_nodes(), coloring, curtime)

        # plt.show()
        # plt.show(block=True)
        fig.savefig(filename)
        # print G.size()
        if return_layout:
            return pos


class TestClassMethods(unittest.TestCase):
    """
    :Purpose:
    	Unittest for testing the methods of the HIV class.
    """

    def setUp(self):
        """
        :Purpose:
            Tests that all models from setup pass inspection.
            ``setUp`` is perfomed before each method.
        """
        self.N_pop = 10000

    def test_NormalAgents(self):
        """Test if all non-IDU,ND,NIDU agents are in the population"""
        print(" ... Test: NormalAgents")
        myNetworkObj = NetworkClass(N=10000, m_0=3)
        for agent in myNetworkObj.NormalAgents:
            agent_sex_type = myNetworkObj.get_agent_characteristic(agent, "Sex Type")
            agent_drug_type = myNetworkObj.get_agent_characteristic(agent, "Drug Type")
            self.assertTrue(
                agent_drug_type not in ["IDU", "NIDU"], "Wrong drug type!%s" % str(agent_drug_type)
            )
            self.assertTrue(agent_sex_type != "MSM", "Wrong sex type!%s" % str(agent_sex_type))

    def test_PartialNetwork(self):
        """Test if all non-IDU,ND,NIDU agents are in the population"""
        print(" ... Test: Network and Agent type consistency")
        myNetworkObj = NetworkClass(N=10000, m_0=3)

        for agent in myNetworkObj.NormalAgents:
            self.assertTrue(agent in myNetworkObj.G.nodes())

        for agent in myNetworkObj.G:
            agent_sex_type = myNetworkObj.get_agent_characteristic(agent, "Sex Type")
            agent_drug_type = myNetworkObj.get_agent_characteristic(agent, "Drug Type")
            self.assertTrue(
                agent_drug_type not in ["IDU", "NIDU"], "Wrong drug type!%s" % str(agent_drug_type)
            )
            self.assertTrue(agent_sex_type != "MSM", "Wrong sex type!%s" % str(agent_sex_type))

    def test_PopulationConsistency(self):
        """Test if Drug users add up"""
        print(" ... Test: Population consistency")
        myNetworkObj = NetworkClass(N=10000, m_0=3)
        CheckSumDrug = (
            len(myNetworkObj.IDU_agents)
            + len(myNetworkObj.NIDU_agents)
            + len(myNetworkObj.ND_agents)
        )
        self.assertTrue(myNetworkObj.PopulationSize == CheckSumDrug)

    def test_HIVConsistency(self):
        """Test HIV consistency"""
        print(" ... Test: Test HIV consistency")
        myNetworkObj = NetworkClass(N=10000, m_0=3)
        for agent in myNetworkObj.Agents:
            HIV_status = myNetworkObj.get_agent_characteristic(agent, "HIV")
            if HIV_status == 1:
                self.assertTrue(agent in myNetworkObj.HIV_agents)

        for agent in myNetworkObj.HIV_agents:
            HIV_status = myNetworkObj.get_agent_characteristic(agent, "HIV")
            self.assertTrue(HIV_status == 1)


if __name__ == "__main__":
    unittest.main()
