# -*- coding: utf-8 -*-
"""
Created on Sat Feb  2 19:46:38 2019

@author: Admin
"""

"""
====================
Parallel Betweenness
====================

Example of parallel implementation of betweenness centrality using the
multiprocessing module from Python Standard Library.

The function betweenness centrality accepts a bunch of nodes and computes
the contribution of those nodes to the betweenness centrality of the whole
network. Here we divide the network in chunks of nodes and we compute their
contribution to the betweenness centrality of the whole network.

This doesn't work in python2.7.13. It does work in 3.6, 3.5, 3.4, and 3.3.

It may be related to this:
https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-multiprocessing-pool-map
"""

#from multiprocessing import set_start_method, Pool
from multiprocessing import Pool
import time
import itertools
import sys
import matplotlib.pyplot as plt
import networkx as nx

import math
import networkx.utils
from networkx.utils import powerlaw_sequence
from operator import itemgetter, attrgetter, methodcaller
from collections import OrderedDict

#jobID = sys.argv[1]
jobID = '04'
noofcores = 8

#set_start_method('forkserver')

def chunks(l, n):
    """Divide a list of nodes `l` in `n` chunks"""
    l_c = iter(l)
    while 1:
        x = tuple(itertools.islice(l_c, n))
        if not x:
            return
        yield x


def _betmap(G_normalized_weight_sources_tuple):
    """Pool for multiprocess only accepts functions with one argument.
    This function uses a tuple as its only argument. We use a named tuple for
    python 3 compatibility, and then unpack it when we send it to
    `betweenness_centrality_source`
    """
    return nx.betweenness_centrality_source(*G_normalized_weight_sources_tuple)


def betweenness_centrality_parallel(G, processes=None):
    """Parallel betweenness centrality  function"""
    p = Pool(noofcores)
    node_divisor = len(p._pool) * 4
    node_chunks = list(chunks(G.nodes(), int(G.order() / node_divisor)))
    num_chunks = len(node_chunks)
    bt_sc = p.map(_betmap,
                  zip([G] * num_chunks,
                      [True] * num_chunks,
                      [None] * num_chunks,
                      node_chunks))
    # Reduce the partial solutions
    bt_c = bt_sc[0]
    for bt in bt_sc[1:]:
        for n in bt:
            bt_c[n] += bt[n]
    p.close()
    p.join()
    return bt_c


def runcent(reffilename,jobID):  
        betweenName = jobID + "betweeness"
        degreeName = jobID + "degree"
        G = nx.read_edgelist("%s_GiantComponent.txt" %reffilename)
        G.remove_edges_from(G.selfloop_edges())
        nx.write_edgelist(G,"%s_processed_GC.txt"%reffilename)
        d = dict(G.degree)

        n = OrderedDict(sorted(d.items(), key=lambda t: t[1], reverse=True))
        with open ('%s.txt'%degreeName, 'w') as fp:
            for p in n.items():
                fp.write("%s:%s\n" % p)
        print("\tParallel version")
        start = time.time()
        bt = betweenness_centrality_parallel(G)
        print("\t\tTime: %.4F" % (time.time() - start))
        
        b_sorted = OrderedDict(sorted(bt.items(), key=lambda t: t[1], reverse=True))
        with open ('%s.txt'%betweenName, 'w') as fp:
            for p in b_sorted.items():
                fp.write("%s:%s\n" % p)


#if __name__ == "__main__":
#    runcent(reffilename,jobID)                
                
                
#        print("\t\tBetweenness centrality for node 0: %.5f" % (bt[0]))
#        print("\tNon-Parallel version")
#        start = time.time()
#        bt = nx.betweenness_centrality(G)
#        print("\t\tTime: %.4F seconds" % (time.time() - start))
#        print("\t\tBetweenness centrality for node 0: %.5f" % (bt[0]))
#    print("")

#    nx.draw(G_ba)
#    plt.show()
