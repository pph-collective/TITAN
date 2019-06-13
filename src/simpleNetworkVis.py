from __future__ import absolute_import
from __future__ import print_function
import matplotlib.pyplot as plt
import networkx as nx
import random as random
import numpy as np
from scipy import stats

num_agents = 10000
graphNum = 0
maxPodSize = 100
G = nx.Graph()
while(num_agents>0):
  graphNum+=1
  lowerLim = min(num_agents+2, maxPodSize)
  curPodN = random.randrange(2,lowerLim)
  num_agents -= curPodN

  print(graphNum, curPodN)
  if (curPodN == 2):
    tempG = nx.random_geometric_graph(curPodN, 1.0)
  else:
    tempG = nx.connected_watts_strogatz_graph(curPodN,2,0.125, tries=200)
  #tempG = nx.random_geometric_graph(curPodN, 0.125)
  while(not nx.is_connected(tempG)):
    tempG = nx.random_geometric_graph(curPodN, 0.125)
  #tempG.graph['ID']=graphNum
  G = nx.disjoint_union(G, tempG)
  
print(nx.number_of_nodes(G))
print(nx.number_connected_components(G))
print(nx.info(G))  
from networkx.algorithms import approximation

components = sorted(nx.connected_component_subgraphs(G), key=len, reverse=True)
bigG = components[0]
outfile = open('networkStats.txt','w')
outfile.write(nx.info(G))

centDict = nx.degree_centrality(G)

outfile.write("\nNumber of connected components: {}\n".format(nx.number_connected_components(G)))
conG = nx.connected_components(G)
tot_nodes = 0
for t in conG:
    thisComp = len(t)
    tot_nodes += thisComp
outfile.write("Average component size: {}\n".format(tot_nodes*1./nx.number_connected_components(G)))
outfile.write("Maximum component size: {}\n".format(nx.number_of_nodes(bigG)))
outfile.write("Degree Histogram: {}\n".format(nx.degree_histogram(G)))
outfile.write("Graph density: {}\n".format(nx.density(G)))
outfile.write("Average node degree centrality: {}\n".format(sum(centDict.values())/len(list(centDict.values()))))
#outfile.write("Average node connectivity: {}\n".format(nx.node_connectivity(G)))
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
#print np.histogram(comps, bins=np.arange(maxPodSize))

print(stats.describe(comps))

import matplotlib.pyplot as plt
plt.hist(comps, bins='auto')  # arguments are passed to np.histogram
plt.title("Histogram with 'auto' bins")
plt.show()





draw = False

if(draw):

  # position is stored as node attribute data for random_geometric_graph
  pos = nx.get_node_attributes(G, 'pos')
  # find node near center (0.5,0.5)
  dmin = 1
  ncenter = 0
  for n in pos:
      x, y = pos[n]
      d = (x - 0.5)**2 + (y - 0.5)**2
      if d < dmin:
          ncenter = n
          dmin = d
  # color by path length from node near center
  p = dict(nx.single_source_shortest_path_length(G, ncenter))

  plt.figure(figsize=(8, 8))
  nx.draw_networkx_edges(G, pos, nodelist=[ncenter], alpha=0.4)
  nx.draw_networkx_nodes(G, pos, nodelist=list(p.keys()),
                         node_size=80,
                         node_color=list(p.values()),
                         cmap=plt.cm.Reds_r)

  plt.xlim(-0.05, 1.05)
  plt.ylim(-0.05, 1.05)
  plt.axis('off')
  plt.show()