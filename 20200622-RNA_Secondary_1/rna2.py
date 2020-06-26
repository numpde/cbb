# RA, 2020-06-22

"""
RNA secondary structure prediction, second algorithm of [2].
Refs: [1, p.209], [3].

[1] Clote & Backofen, Computational Molecular Biology: An Introduction, 2000, Wiley
[2] Nussinov et al., Algorithms for loop matching, SIAM, 1978
[3] Nussinov and Jacobson, Fast algorithm for [..] RNA, PNAS, 1980
[4] https://en.wikipedia.org/wiki/Potato_spindle_tuber_viroid
"""

import os
import io
import numpy as np
import pandas as pd
import networkx as nx

from networkx.drawing.nx_pydot import graphviz_layout

from inclusive import range
from itertools import product, cycle, count
from collections import defaultdict
from tcga.complements import dna_to_dna
from tcga.utils import download
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


# Pairing energy [1, p.209]
def a(x, y):
    xy = x + y
    if (xy == "GC") or (xy == "CG"):
        return -5
    if (xy == "AT") or (xy == "TA"):
        return -4
    if (xy == "GT") or (xy == "TG"):
        return -1
    return np.inf


def pair_up1(s):
    def min_e(s, memory={}):
        """
        Minimize total pairing energy.
        Refs: [1, p.209], [3, Second algorithm].
        """

        if (len(s) <= 1):
            return (0, "." * len(s))

        if s not in memory:
            # Energy minimum
            m = 0
            cand = "." * len(s)

            if (len(s) >= 3):
                e = a(s[0], s[-1])
                if (e < np.inf):
                    (c, p) = min_e(s[1:-1])
                    if (c + e <= m):
                        m = c + e
                        cand = "(" + p + ")"

            for j in range(1, len(s)):
                (c1, p1) = min_e(s[:j])
                (c2, p2) = min_e(s[j:])
                if (c1 + c2 <= m):
                    m = c1 + c2
                    cand = p1 + p2

            memory[s] = (m, cand)

        return memory[s]

    return min_e(s)


def pair_up2(s):
    # Pairwise energies
    E = {(i, j): a(x, y) for (i, x) in enumerate(s) for (j, y) in enumerate(s)}
    M = defaultdict(int)
    S = defaultdict(str)
    K = defaultdict(list)
    for i in range(len(s)):
        S[(i, i)] = "."
    for n in range(1, len(s)):
        for (i, j) in zip(count(0), range(n, len(s))):
            assert (0 <= i < (i + n) == j < len(s))
            m = 0
            if (abs(i - j) >= 2):
                e = E[(i, j)] + M[(i + 1, j - 1)]
                if (e <= m):
                    m = e
                    M[(i, j)] = m
                    S[(i, j)] = "(" + S[(i + 1, j - 1)] + ")"
                    K[(i, j)] = [(i + 1, j - 1)]
            for k in range(i, j):
                assert (i <= k < j)
                [p, q] = [(i, k), (k + 1, j)]
                e = M[p] + M[q]
                if (e <= m):
                    m = e
                    M[(i, j)] = m
                    S[(i, j)] = S[p] + S[q]
                    K[(i, j)] = [p, q]

    def brackets(ij):
        if ij not in K:
            assert ij[0] == ij[1]
            return "."
        if (len(K[ij]) == 1):
            return "(" + brackets(K[ij][0]) + ")"
        if (len(K[ij]) == 2):
            return brackets(K[ij][0]) + brackets(K[ij][1])
        raise RuntimeError

    ij = (0, len(s) - 1)
    B = brackets(ij)
    assert (S[ij] == B)

    return (M[ij], B)


def construct_graph(S: str, C: str) -> nx.Graph:
    """
    S is the nucleotide string.
    C is the bracket pattern.
    """
    g = nx.Graph()
    stack = []
    for (i, s) in enumerate(S):
        g.add_node(i, n=s)
    for (i, j) in zip(range(len(S)), list(range(1, len(S))) + [0]):
        g.add_edge(i, j, type="bb")
    for (i, (s, c)) in enumerate(zip(S, C)):
        if (c == "("):
            stack += [(i, s)]
        if (c == ")"):
            (j, _) = stack.pop()
            g.add_edge(i, j, type="bp")
    return g


def view_graph(g):
    with Plox() as px:
        bp = [(a, b) for (a, b, d) in g.edges.data(data='type') if (d == 'bp')]
        bb = [(a, b) for (a, b, d) in g.edges.data(data='type') if (d == 'bb')]

        pos = graphviz_layout(g, prog="sfdp")
        # pos = nx.planar_layout(g)
        pos = nx.spring_layout(g, k=10, pos=pos, iterations=100)

        params = dict(node_size=2)
        nx.draw_networkx_nodes(g, pos=pos, **params)
        nx.draw_networkx_edges(g, pos=pos, ax=px.a, edgelist=bp, edge_color='r', **params)
        nx.draw_networkx_edges(g, pos=pos, ax=px.a, edgelist=bb, edge_color='k', **params)

        PARAM['out_fig'].parent.mkdir(parents=True, exist_ok=True)
        px.f.savefig(PARAM['out_fig'])


def main():
    viroid_fasta = download(PARAM['viroid']).to(rel_path="cache/download").now.text
    pstg = SeqIO.read(io.StringIO(viroid_fasta), format='fasta')
    S = pstg.seq[0:]
    # (e, C) = pair_up1(S)
    # print(e)
    # print(S)
    # print(C)
    (e, C) = pair_up2(S)
    print(e)
    print(S)
    print(C)
    # g = construct_graph(S, C)
    # view_graph(g)


def test():
    pass


if __name__ == '__main__':
    os.nice(19)
    test()
    main()
