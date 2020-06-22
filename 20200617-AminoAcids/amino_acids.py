# RA, 2020-06-17

"""
Display amino acids in a tree, colored by property.

(Unfinished)
"""

import numpy as np
import pandas as pd
import networkx as nx

from plox import Plox
from pathlib import Path
from itertools import groupby
from contextlib import contextmanager
from tcga.codons.tables import standard as rna_codons
from more_itertools import last
from networkx.drawing.nx_pydot import graphviz_layout

ROOT = Path(__file__).parent
mkdir = (lambda x: (x.mkdir(exist_ok=True, parents=True)) or x)

PARAM = {
    'AA': ROOT / "aa_prop1.csv",

    'out_dir': ROOT / "figs",
}


def hierarchical_codons():
    return {
        n1: {
            n2: {
                n3: rna_codons[n1 + n2 + n3]
                for (n3, g3) in groupby(sorted(list(g2)), key=(lambda x: x[2]))
            }
            for (n2, g2) in groupby(sorted(list(g1)), key=(lambda x: x[1]))
        }
        for (n1, g1) in groupby(sorted(rna_codons), key=(lambda x: x[0]))
    }


# def grouped_codons():
#     {'a': {'a': {'ag': 'K', 'cu': 'N'}, 'c': {'acgu': 'T'},
#            'g': {'ag': 'R', 'cu': 'S'}, 'u': {'acu': 'I', 'g': 'M'}},
#      'c': {'a': {'ag': 'Q', 'cu': 'H'}, 'c': {'acgu': 'P'},
#            'g': {'acgu': 'R'}, 'u': {'acgu': 'L'}},
#      'g': {'a': {'ag': 'E', 'cu': 'D'}, 'c': {'acgu': 'A'},
#            'g': {'acgu': 'G'}, 'u': {'acgu': 'V'}},
#      'u': {'a': {'ag': 'X', 'cu': 'Y'}, 'c': {'acgu': 'S'},
#            'g': {'a': 'X', 'cu': 'C', 'g': 'W'}, 'u': {'ag': 'L', 'cu': 'F'}}}


def get_aa():
    return pd.read_csv(PARAM['AA'], sep='\t', comment='#')


def make_aa_graph():
    g = nx.DiGraph()
    o = "___"
    g.add_node(o, kind="origin")

    C = {
        F"{n1}__": {
            F"{n1}{n2}_": {
                F"{n1}{n2}{n3}": list(g3)
                for (n3, g3) in g2.items()
            }
            for (n2, g2) in g1.items()
        }
        for (n1, g1) in hierarchical_codons().items()
    }

    for (n1, g1) in C.items():
        g.add_edge(o, n1)
        for (n2, g2) in g1.items():
            g.add_edge(n1, n2)
            for (n3, g3) in g2.items():
                g.add_edge(n2, n3)

                aa = rna_codons[n3]
                g.add_node(aa, kind="aa")
                g.add_edge(n3, aa)

    return g


@contextmanager
def visualize_graph(g: nx.DiGraph) -> Plox:
    with Plox() as px:
        nodes_0 = ["___"]
        nodes_1 = [n for (n, k) in g.nodes(data='kind') if (k != "aa") and (n.count('_') == 2)]
        nodes_2 = [n for (n, k) in g.nodes(data='kind') if (k != "aa") and (n.count('_') == 1)]
        nodes_3 = [n for (n, k) in g.nodes(data='kind') if (k != "aa") and (n.count('_') == 0)]
        nodes_aa = [n for (n, k) in g.nodes(data='kind') if (k == "aa")]

        pos = graphviz_layout(g, prog="twopi", root='___')
        # pos = nx.spring_layout(g, pos=pos)
        # pos = nx.shell_layout(g, nlist=[nodes_0, nodes_1, nodes_2, nodes_3, nodes_aa])
        # pos = nx.planar_layout(g)
        # pos = nx.kamada_kawai_layout(g)
        # pos = nx.spring_layout(g, pos=pos, k=10, iterations=10, threshold=1e-8)

        nx.draw_networkx_edges(g, pos=pos)

        nx.draw_networkx_nodes(g, pos=pos, nodelist=(nodes_0 + nodes_1 + nodes_2 + nodes_3))
        nx.draw_networkx_nodes(g, pos=pos, nodelist=nodes_aa, node_color='r')

        labels = {n: {'aa': n, 'origin': "-", None: last(n.strip('_'), None)}[k] for (n, k) in g.nodes(data='kind')}
        nx.draw_networkx_labels(g, pos=pos, labels=labels)

        yield px


def main():
    g = make_aa_graph()

    with visualize_graph(g) as px:
        px.f.savefig(mkdir(PARAM['out_dir']) / "aa.png")

    # aa = get_aa()


if __name__ == '__main__':
    main()
