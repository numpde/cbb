# RA, 2020-06-22

"""
RNA secondary structure prediction, second algorithm of [2].
Refs: [1, p.209], [3].

See [5] for thermodynamic considerations, in particular
"stability numbers" for pairs, hairpin loops, interior loops and bulges,
relative to the single strand state in units of A-U base pair energy.
Note on [5, p.364]:
    "...the practice of simply maximizing base pairs to predict
    secondary structure must be tempered by also minimizing
    the number of loops."

[1] Clote & Backofen, Computational Molecular Biology: An Introduction, 2000, Wiley
[2] Nussinov et al., Algorithms for loop matching, SIAM, 1978
[3] Nussinov and Jacobson, Fast algorithm for [..] RNA, PNAS, 1980
[4] https://en.wikipedia.org/wiki/Potato_spindle_tuber_viroid
[5] Tinoco, Uhlenbeck and Levine, Estimation of secondary structure of RNA, Nature, 1971
"""

import os
import io
import numpy as np
import pandas as pd
import networkx as nx

from networkx.drawing.nx_pydot import graphviz_layout

from inclusive import range
from itertools import product, cycle, count, repeat
from collections import defaultdict
from tcga.complements import dna_to_dna, dna_to_rna
from tcga.utils import download, First
from Bio import SeqIO
from plox import Plox
from pathlib import Path
from diskcache import Cache

ROOT_PATH = Path(__file__).parent

PARAM = {
    'viroid': "https://www.ncbi.nlm.nih.gov/search/api/download-sequence/?db=nuccore&id=V01465.1",
    'out_fig': ROOT_PATH / "figs/rna2.png",
}

diskcache = Cache(ROOT_PATH / "cache/compute/UV/").memoize(expire=None)


# Pairing energy, modified from [1, p.209]
def a(x, y):
    xy = x + y
    if (xy == "GC") or (xy == "CG"):
        return -6
    if (xy == "AU") or (xy == "UA"):
        return -4
    if (xy == "GU") or (xy == "UG"):
        return -1
    return np.inf

class R:
    def __init__(self, e=np.inf, s="-"):
        (self.e, self.s) = (e, s)

    def __add__(self, other):
        return R(self.e + other.e, self.s + other.s)

    def __lt__(self, other):
        return (self.e < other.e)

    def as_tuple(self):
        return (self.e, self.s)


def pair_up2(s) -> R:
    # Pairwise energies
    E = {
        (i, j): a(x, y)
        for (i, x) in enumerate(s) for (j, y) in enumerate(s)
        if (6 <= min(abs(i - j), abs(i - (j - len(s))), abs((i - len(s)) - j)))
    }

    mem = defaultdict(R)

    def rec(i, j) -> R:
        if (i == j):
            return R(0, ".")
        if (i, j) in mem:
            return mem[(i, j)]
        if (i, j) in E:
            mem[(i, j)] = R(E[(i, j)], "(") + rec(i + 1, j - 1) + R(0, ")")
        for (p, q) in zip(zip(repeat(i), range(i, j)), zip(count(i + 1), repeat(j))):
            mem[(i, j)] = min(mem[(i, j)], rec(*p) + rec(*q))
        return mem[(i, j)]

    return rec(0, len(s) - 1)


def pairs(B: str):
    """
    For a string like B = '.(().)' yield the pairs (i, j)
    of indices of matching brackets, i.e. (1, 5), (2, 3).
    """
    stack = []
    for (i, b) in enumerate(B):
        if (b == "("):
            stack.append(i)
        if (b == ")"):
            yield (stack.pop(), i)


def construct_graph(S: str, B: str) -> nx.Graph:
    """
    S is the nucleotide string.
    B is the bracket pattern.
    """
    g = nx.Graph()
    stack = []
    for (i, s) in enumerate(S):
        g.add_node(i, n=s)
    for (i, j) in zip(range(len(S)), list(range(1, len(S))) + [0]):
        g.add_edge(i, j, type="bb")
    for (i, j) in pairs(B):
        g.add_edge(i, j, type="bp")
    return g


def view_graph(g):
    with Plox() as px:
        params = dict(node_size=2)

        bp = [(a, b) for (a, b, d) in g.edges.data(data='type') if (d == 'bp')]
        bb = [(a, b) for (a, b, d) in g.edges.data(data='type') if (d == 'bb')]

        for (a, b, d) in g.edges.data(data='type'):
            g.edges[(a, b)]['weight'] = {'bp': 1, 'bb': 10}[d]

        pos = graphviz_layout(g, prog="sfdp")
        pos = nx.spring_layout(g, k=10, pos=pos, iterations=1000, threshold=1e-8, weight='weight')

        # pos = graphviz_layout(g, prog="circo")
        # for i in range(1000000):
        #     for (a, b, d) in g.edges(data='type'):
        #         k = {'bp': 50, 'bb': 10}[d]
        #         va = np.asarray(pos[a])
        #         vb = np.asarray(pos[b])
        #         l = np.linalg.norm(va - vb)
        #         f = 0.9 if (l > k) else 1.01
        #         (va, vb) = (vb + (va - vb) * f, va + (vb - va) * f)
        #         pos[a] = tuple(va)
        #         pos[b] = tuple(vb)
        #
        #     if not (i % 1000):
        #         import matplotlib.pyplot as plt
        #
        #         px.f.clear()
        #         nx.draw_networkx_nodes(g, pos=pos, **params)
        #         nx.draw_networkx_edges(g, pos=pos, ax=px.a, edgelist=bp, edge_color='b', **params)
        #         nx.draw_networkx_edges(g, pos=pos, ax=px.a, edgelist=bb, edge_color='k', **params)
        #         # nx.draw_n
        #
        #         plt.ion()
        #         plt.show()
        #         plt.pause(0.1)
        #
        #
        # exit()


        nx.draw_networkx_nodes(g, pos=pos, **params)
        nx.draw_networkx_edges(g, pos=pos, ax=px.a, edgelist=bp, edge_color='b', **params)
        nx.draw_networkx_edges(g, pos=pos, ax=px.a, edgelist=bb, edge_color='k', **params)

        nx.draw_networkx_nodes(g, pos=pos, ax=px.a, nodelist=[min(g.nodes)], node_color='g', **params)
        nx.draw_networkx_nodes(g, pos=pos, ax=px.a, nodelist=[max(g.nodes)], node_color='r', **params)

        PARAM['out_fig'].parent.mkdir(parents=True, exist_ok=True)
        px.f.savefig(PARAM['out_fig'])


def main():
    viroid_fasta = download(PARAM['viroid']).to(rel_path="cache/download").now.text
    pstg = SeqIO.read(io.StringIO(viroid_fasta), format='fasta')

    S = pstg.seq[0:]
    S = First(dna_to_dna).then(dna_to_rna)(S)

    (e, B) = pair_up2(S).as_tuple()
    print("Folding energy:", e)
    print("Primary:  ", S)
    print("Secondary:", B)

    print(
        "Number of base pairs:",
        {
            p: sum(
                (S[i] + S[j]) in p
                for (i, j) in pairs(B)
            )
            for p in ["GC/CG", "AU/UA", "GU/UG"]
        }
    )

    g = construct_graph(S, B)
    view_graph(g)


def test():
    pass


if __name__ == '__main__':
    os.nice(19)
    test()
    main()
