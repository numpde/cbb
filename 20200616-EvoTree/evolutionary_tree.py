# RA, 2020-06-16

"""
Proof-of-concept reconstruction of a tree from distances.

May require
    sudo apt install graphviz

Follows
    Thm 4.5 in [1, p.154] and
    Unweighted Pair Group Method with Arithmetic mean -- Alg 4.1 in [1, on p.149].

[1] Clote & Backofen, Computational Molecular Biology: An Introduction, 2000, Wiley
[2] Farris, Estimating phylogenetic trees from distance matrices, The American Naturalist, 1972
[3] KKBM, Calculation of evolutionary trees from sequence data, PNAS, 1979
"""

import ast
from typing import Any, Tuple

import numpy as np
import pandas as pd
import networkx as nx

from plox import Plox
from pathlib import Path
from itertools import product
from networkx.drawing.nx_pydot import graphviz_layout

ROOT = Path(__file__).parent
mkdir = (lambda x: (x.mkdir(exist_ok=True, parents=True)) or x)

PARAM = {
    'out_path': mkdir(ROOT / "figs"),
}


def random_tree() -> nx.DiGraph:
    g: nx.DiGraph
    rs = np.random.RandomState(seed=1)
    g = nx.bfs_tree(nx.random_tree(n=16, seed=1), source=0)
    # nx.set_edge_attributes(g, 1, name='len')
    nx.set_edge_attributes(g, dict(zip(g.edges, rs.uniform(1, 3, size=len(g.edges)))), name='len')
    return g


def save_tree(g: nx.DiGraph, filename: Path):
    with Plox() as px:
        nx.draw(g, pos=graphviz_layout(g, prog="dot"), ax=px.a, with_labels=True)
        px.f.savefig(filename)
    return g


def farris(d: pd.DataFrame, r, leaves):
    # Average distance between r and the leaves
    c = d[r][leaves].mean()
    # Transformed distances according to Farris / KKBM
    e = d.copy(deep=True).loc[leaves, leaves]
    for (i, j) in product(e.index, e.columns):
        # This is Eq (7) in [3]
        e.loc[i, j] = ((d.loc[i, j] - d.loc[i, r] - d.loc[r, j]) + 2 * c) if (i != j) else 0
    return e


def upgma(d: pd.DataFrame):
    d = d.copy(deep=True)

    for i in d.index:
        d.loc[i, i] = np.nan

    for n in range(len(d) - 1):
        (i, j) = np.nonzero((d == d.min().min()).to_numpy())
        (i, j) = next(zip(d.index[i], d.columns[j]))
        k = "({}, {})".format(i, j)
        r = pd.DataFrame({k: (d[i] + d[j]) / 2})
        d = pd.concat([d, r], axis=1)
        d = pd.concat([d, r.T], axis=0)
        d = d.drop(columns=[i, j], index=[i, j])

    return tuple(ast.literal_eval(k))


def uncluster(d: pd.DataFrame, cc) -> Tuple[Any, nx.DiGraph]:
    g: nx.DiGraph
    g = nx.DiGraph()
    g.add_node(cc)

    if type(cc) is not tuple:
        return (cc, g)

    assert (len(cc) == 2)

    ((u, U), (v, V)) = tuple(uncluster(d, c) for c in cc)

    g = nx.union_all([g, U, V])
    g.add_edge(cc, u)
    g.add_edge(cc, v)

    return (cc, g)


def get_leaves(g: nx.DiGraph):
    return [x for x in g.nodes() if g.in_degree(x) and not g.out_degree(x)]


def main():
    g = random_tree()
    leaves = get_leaves(g)
    save_tree(g, PARAM['out_path'] / "original_tree.png")

    d = pd.DataFrame(dict(nx.all_pairs_dijkstra_path_length(g.to_undirected(), weight='len')))

    e = upgma(farris(d, r=0, leaves=leaves))

    (r, h) = uncluster(d, e)

    save_tree(h, PARAM['out_path'] / "restored_tree.png")


if __name__ == '__main__':
    main()
