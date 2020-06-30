# RA, 2020-06-29

"""
A networkx based reimplementation of rna1.py.
"""

import os

import numpy as np
import networkx as nx

from rna2 import get_pstg_seq, a


def pair_up(s):
    g = nx.DiGraph()
    g.add_node((0, 0))
    for j in range(len(s)):
        for (i, n) in list(g.nodes):
            g.add_edge((i, n), (j, n + 1), weight=0)
            if (n > 0):
                g.add_edge((i, n), (j, n - 1), weight=a(s[i], s[j]))
    for (i, n) in list(g.nodes):
        if (n == 0):
            g.add_edge((i, n), (len(s), 0), weight=0)
    print("Graph OK")
    return nx.bellman_ford_path_length(g, (0, 0), (len(s), 0))


def main():
    s = get_pstg_seq()[0:100]
    print(pair_up(s))


if __name__ == '__main__':
    os.nice(19)
    main()
