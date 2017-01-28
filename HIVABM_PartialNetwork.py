#!/usr/bin/env python2.3
# -*- coding: utf-8 -*-
"""
Module that carries all definitions 
for the ABM of drug use and HIV. \n

Author:		Lars Seemann \n
Email:		lseemann@uh.edu \n 
Date:		2011-01-28 \n

Copyright (c) 2010, under the Simplified BSD License. \n
For more information on FreeBSD see: http://www.opensource.org/licenses/bsd-license.php \n
All rights reserved.
"""
__author__="Lars Seemann (lseemann@uh.edu)"

import os
import random
import copy
import itertools
from copy import deepcopy
import unittest
import numpy as np

try: from HIVABM_Population import PopulationClass
except ImportError:
    raise ImportError("Can't import PopulationClass")

try: import matplotlib.pyplot as plt
except ImportError:
    raise ImportError("You must install matplotlib (http://matplotlib.sourceforge.net/)")
#try: import networkx as nx
#except ImportError:
#    raise ImportError("You must install NetworkX (http://networkx.lanl.gov/) ")

#def save_adjlist(graph, dir_prefix, time):
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
    #nx.write_adjlist(graph,(dir_prefix +'/'+'test_adjlist_at_time_%d.txt'%time))
    if not os.path.isdir(dir_prefix): os.mkdir(dir_prefix)
    OutFileName = (dir_prefix + '/' +'Interactions_at_time_%d.txt'%time)
    outfile = open(OutFileName,'w')
    outfile.write('Agent 1\ttime\tAgent 2\n')
    not_singletons = []
    
    for e in graph.edges_iter():
        not_singletons.append(e[0])
        not_singletons.append(e[1])
        outfile.write('%d\t%d\t%d\n'%(e[0],time,e[1]))
    singletons = list(set(range(N_pop)).difference(set(not_singletons)))
    for singleton in singletons:
        outfile.write('%d\t%d\n'%(singleton, time))
        
def _random_subset(seq,m):
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
    G.name = "my_erdos_renyi_binomial_random_graph(%s,%s)"%(n,p)
    if p <= 0:
        return G
    if p >= 1:
        return complete_graph(n,create_using = G)

    if not seed is None:
        random.seed(seed)

    if G.is_directed():
        edges = itertools.permutations(node_list,2)
    else:
        edges=itertools.combinations(node_list,2)

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

    if m < 1 or  m >=n:
        raise nx.NetworkXError(\
            "Barabási-Albert network must have m>=1 and m<n, m=%d,n=%d"%(m,n))

    if node_list is not None:
        if n != len(node_list):
            raise nx.NetworkXError(\
                "node_list must have smae length as n! len(node_list)=%d,n=%d"%(len(node_list),n))
        else:
            random.shuffle(node_list)              # Shuffle node_list
            targets = node_list[:m]
    else:
        targets=random.shuffle(list(range(m)))     # Target nodes for new edges

    if seed is not None:
        random.seed(seed)   

    G = nx.Graph()                                 # networkX graph
    G.add_nodes_from(targets[:m])                  # Add m initial nodes (m0 in barabasi-speak)  
    G.name = "mod_barabasi_albert_graph(%s,%s)"%(n,m)
    repeated_nodes = []                            # List of existing nodes, with nodes repeated once for each adjacent edge 
    Num_source = m                                 # Start adding the other n-m nodes. The first node is m.
    while Num_source < n: 
        source = node_list[Num_source]
        G.add_edges_from(zip([source]*m,targets)) # Add edges to m nodes from the source.
        repeated_nodes.extend(targets)            # Add one node to the list for each new edge just created.
        repeated_nodes.extend([source]*m)         # And the new node "source" has m edges to add to the list.
        # Now choose m unique nodes from the existing nodes 
        # Pick uniformly from repeated_nodes (preferential attachement) 
        targets = _random_subset(repeated_nodes,m)
        Num_source += 1
    return G

class NetworkClass(PopulationClass):

    def __init__(self, N = 10000, m_0 = 1, network_type='scale_free'):
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
        if type(N) is not int:			
            raise ValueError(('Population size must be integer,\
                                      n = %s, not %s')%(string(N), type(N)))	
        else: pass
        if m_0 not in range(10):
            raise ValueError('m_0 must be integer smaller than 10')
        else: self.m_0 = m_0

        PopulationClass.__init__(self, n = N)	# Create population
        self.NormalAgents = []
        for agent in self.Agents:
            if (agent not in self.IDU_agents and 
                agent not in self.NIDU_agents and 
                agent not in self.MSM_agents):
                self.NormalAgents.append(agent)
        self.NetworkSize = len(self.NormalAgents)
        if network_type=='scale_free':
            # scale free Albert Barabsai Graph
            self.G = my_barabasi_albert_graph(self.NetworkSize,m_0, 
                                              node_list = self.NormalAgents)
        elif network_type=='binomial':
            self.G = my_erdos_renyi_binomial_random_graph(
                node_list=self.NormalAgents, 
                p=0.5, 
                seed=None, 
                directed=False)
        else:
            raise ValueError("Invalid network type! %s"%str(network_type))


    def plot_DegreeDistribution(self, graph):
        """ 
        Plot the node degree distribution of the graph. \n 
        INPUT: networkX graph
        """
        print("\tNetwork size = "+str(self.NetworkSize))
        degreeList = np.array( graph.degree().values() )	
        x_degree = np.arange(max(degreeList)+1)	# include 0 and max
        degreehist = np.bincount(degreeList)
        ix = degreehist != 0			# delete zeros
        degreehist = degreehist[ix]
        x_degree = x_degree[ix]

        # print node degree distribution	
        plt.loglog(x_degree, degreehist, 'bo')
        plt.title('Node degree distribution')
        plt.ylabel('Frequency')
        plt.xlabel('Node degree')
        plt.axis('equal')
        plt.axis([min(degreeList)-(0.05*max(degreeList)), 
                  max(degreeList)+(0.05*max(degreeList)),
                  0.9, max(degreehist)+(0.05* max(degreehist))])
        plt.show()

    def get_AdjacencyList(self):
        """ Return the adjacency list of the graph. """
        return self.G.adjacency_list()

    def get_Graph(self):
        """
        Return random assortative graph produced by ``set_assortative_graph``.
        """		
        return self.G

    def visualize_network(self, graph, coloring='Sex Type', pos=None, 
                          return_layout=0, node_size=None):
        """
        :Purpose:
            Visualize the network using the spring layout (default). \n

        :Input: 
            graph : networkX graph
        """
        print("Plotting...")
        plt.figure(figsize=(8,8))
        # with nodes sized by degree
        # node_color=[float(H.degree(v)) for v in H]
        # layout:
        if not pos:
            pos=nx.spring_layout(graph, iterations=30)
            #pos = nx.circular_layout(graph)
            #pos=nx.shell_layout(graph)
            #pos=nx.random_layout(graph)
            #pos=nx.spectral_layout(graph)

        edge_color = 'k'
        node_shape = 'o'
        # node color to indicate sex type
        node_color = []
        if coloring=='Sex Type':
            for v in graph:
                tmp_sextype = self.get_agent_characteristic(v, 'Sex Type')
                if tmp_sextype == 'HM':
                    node_color.append('b')
                elif tmp_sextype == 'HF':
                    node_color.append('g')
                elif tmp_sextype == 'WSW':
                    node_color.append('c')
                elif tmp_sextype == 'MSM':
                    node_color.append('r')
                else:
                    raise ValueError("Check agents %s drug %s"%(v, tmp_sextype))
        elif coloring == 'Drug Type':
            for v in graph:
                tmp_drugtype = self.get_agent_characteristic(v, 'Drug Type')
                if tmp_drugtype == 'ND':
                    node_color.append('g')
                elif tmp_drugtype == 'NIDU':
                    node_color.append('b')
                elif tmp_drugtype == 'IDU':
                    node_color.append('r')
                else:
                    raise ValueError("Check agents %s drug type %s"%(v, tmp_drugtype))
        elif coloring == 'AIDS HIV':
            for v in graph:
                tmp_aids = self.get_agent_characteristic(v, 'AIDS')
                tmp_hiv = self.get_agent_characteristic(v, 'HIV')
                if tmp_hiv == 1:
                    node_color.append('r')
                elif tmp_aids == 1:
                    node_color.append('b')
                else:
                    node_color.append('g')
        else:
            raise ValueError("coloring value invalid!\n%s\n \
            Only 'Sex Type','Drug Type', and 'AIDS HIV' allowed!"%str(coloring))

        # node size indicating node degree
        NodeSize=[]
        if node_size:
            for v in graph:
                NodeSize.append(100)
        else:
            for v in graph:
                NodeSize.append((10*graph.degree(int(v)))**(1.0))


        # draw:
        nx.draw(graph, pos, 
                node_size = NodeSize, 
                node_color = node_color,
                node_shape = node_shape,
                edge_color = edge_color,
                alpha=0.4, with_labels=False,
                width=3.0)
        
        #nx.draw_networkx_nodes(self.graph,pos,node_size=NodeSize)
        #nx.draw_networkx_edges(self.graph,pos,alpha=0.4)
        
        plt.axis('equal')
        plt.axis('off')
        filename="Network_%d_%s.png"%(graph.size(),coloring)
        plt.savefig(filename)
        plt.show()
        
        print graph.size()
        if return_layout:
            return pos
        
def main():
    myNetworkObj = ScaleFreeNetworkClass(N=1000, m_0 = 1)
    G = myNetworkObj.get_Graph()
    myNetworkObj.visualize_network(G)
    myNetworkObj.plot_DegreeDistribution(G)

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
        print " ... Test: NormalAgents"
        myNetworkObj = NetworkClass(N=10000, m_0 = 3)
        for agent in myNetworkObj.NormalAgents:
            agent_sex_type = myNetworkObj.get_agent_characteristic(agent,'Sex Type')
            agent_drug_type = myNetworkObj.get_agent_characteristic(agent,'Drug Type')
            self.assertTrue(agent_drug_type not in ['IDU','NIDU'], "Wrong drug type!%s"%str(agent_drug_type))
            self.assertTrue(agent_sex_type != 'MSM',"Wrong sex type!%s"%str(agent_sex_type))

    def test_PartialNetwork(self):
        """Test if all non-IDU,ND,NIDU agents are in the population"""
        print " ... Test: Network and Agent type consistency"
        myNetworkObj = NetworkClass(N=10000, m_0 = 3)

        for agent in myNetworkObj.NormalAgents:
            self.assertTrue(agent in myNetworkObj.G.nodes())

        for agent in myNetworkObj.G:
            agent_sex_type = myNetworkObj.get_agent_characteristic(agent,'Sex Type')
            agent_drug_type = myNetworkObj.get_agent_characteristic(agent,'Drug Type')
            self.assertTrue(agent_drug_type not in ['IDU','NIDU'], "Wrong drug type!%s"%str(agent_drug_type))
            self.assertTrue(agent_sex_type != 'MSM',"Wrong sex type!%s"%str(agent_sex_type))

    def test_PopulationConsistency(self):
        """Test if Drug users add up"""
        print " ... Test: Population consistency"
        myNetworkObj = NetworkClass(N=10000, m_0 = 3)
        CheckSumDrug = len(myNetworkObj.IDU_agents)+len(myNetworkObj.NIDU_agents)+len(myNetworkObj.ND_agents)
        self.assertTrue(myNetworkObj.PopulationSize == CheckSumDrug)

    def test_HIVConsistency(self):
        """Test HIV consistency"""
        print " ... Test: Test HIV consistency"
        myNetworkObj = NetworkClass(N=10000, m_0 = 3)
        for agent in myNetworkObj.Agents:
            HIV_status = myNetworkObj.get_agent_characteristic(agent,'HIV')
            if HIV_status == 1:
                self.assertTrue(agent in myNetworkObj.HIV_agents)

        for agent in myNetworkObj.HIV_agents:
            HIV_status = myNetworkObj.get_agent_characteristic(agent,'HIV')
            self.assertTrue(HIV_status == 1)


if __name__=='__main__':
    unittest.main()
