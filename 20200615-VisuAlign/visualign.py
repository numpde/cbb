# RA, 2020-06-15

"""
Make a matrix plot of one sequence against another.

[1] A human protein related to yeast Cdc6p
    R. S. Williams, R. V. Shohet, and B. Stillman
    "Conservation of structure among proteins involved in initiation suggests
    that fundamental features of replication complexes are maintained in all eukaryotes.
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC19260/
"""

import io
import re
import numpy as np
import pandas as pd
import requests  # pip install requests[security]
import matplotlib.pyplot as plt

from typing import Iterator
from pathlib import Path
from functools import lru_cache as ramcache
from diskcache import Cache
from itertools import product

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

ROOT_PATH = Path(__file__).parent
diskcache = Cache(ROOT_PATH / "cache/UV/").memoize(expire=None)

PARAM = {
    'proteins': {
        # # https://www.uniprot.org/uniprot/P12866
        # STE6, Alpha-factor-transporting ATPase, Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (Baker's yeast)
        # Reviewed, Annotation score: 5/5, Experimental evidence at protein level
        # 'STE6 yeast': "https://www.uniprot.org/uniref/UniRef90_P12866.fasta",

        # https://www.uniprot.org/uniprot/P09119
        # CDC6, Cell division control protein 6, Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (Baker's yeast)
        # Reviewed, Annotation score: 5/5, Experimental evidence at protein level
        'CDC6 yeast': "https://www.uniprot.org/uniref/UniRef90_P09119.fasta",

        # # https://www.uniprot.org/uniprot/A0A024R1S2
        # # CDC6, Cell division control protein (unreviewed)
        # Unreviewed, Annotation score: 2/5, Protein inferred from homology
        # 'A0A024R1S2_HUMAN': "https://www.uniprot.org/uniprot/A0A024R1S2.fasta",

        # https://www.uniprot.org/uniprot/Q99741
        # CDC6, Cell division control protein 6 homolog, Homo sapiens (Human)
        # Reviewed, Annotation score: 5/5, Experimental evidence at protein level
        'CDC6 human': "https://www.uniprot.org/uniref/UniRef90_Q99741.fasta",

        # https://www.uniprot.org/uniprot/O89033
        # Cdc6, Cell division control protein 6 homolog, Mus musculus (Mouse)
        # Reviewed, Annotation score: 5/5, Experimental evidence at protein level
        'CDC6 mouse': "https://www.uniprot.org/uniref/UniRef90_O89033.fasta",

        # https://www.uniprot.org/uniprot/F1R1N7
        # cdc6, Cell division control protein, Danio rerio (Zebrafish) (Brachydanio rerio)
        # Unreviewed, Annotation score: 3/5, Experimental evidence at protein level
        'CDC6 zebrafish': "https://www.uniprot.org/uniprot/F1R1N7.fasta",

        # https://www.uniprot.org/uniprot/Q9VSM9
        # Cdc6, Cell division control protein, Drosophila melanogaster (Fruit fly)
        # Unreviewed, Annotation score: 3/5, Experimental evidence at protein level
        'CDC6 fruit fly': "https://www.uniprot.org/uniprot/Q9VSM9.fasta",

        # https://www.uniprot.org/uniprot/O82387
        # CDC6, Cell division control protein 6 homolog, Arabidopsis thaliana (Mouse-ear cress)
        # Reviewed, Annotation score: 3/5, Experimental evidence at transcript leveli
        'CDC6 arabidopsis': "https://www.uniprot.org/uniprot/O82387.fasta",
    },

    'blosum': ROOT_PATH / "../20200608-Downloads/blosum/blosum80.mat",

    'figdest': ROOT_PATH / "figs",
    'savefig': dict(bbox_inches='tight', pad_inches=0.05, transparent=False, dpi=30),
}


# expire: seconds until arguments expire (default None, no expiry)
@diskcache
def download(url):
    print(F"Downloading {url}")
    return requests.get(url).text


@ramcache
def read_blosum() -> pd.DataFrame:
    with PARAM['blosum'].open('r') as fd:
        s = '\n'.join(re.sub(r"\s+", ' ', line).rstrip() for line in fd.readlines() if not line.strip().startswith('#'))
    blosum = pd.read_csv(io.StringIO(s), sep=' ', index_col=0)
    assert np.array_equal(blosum, blosum.T)
    return blosum


def get_seq() -> Iterator[SeqRecord]:
    for (name, url) in PARAM['proteins'].items():
        for seq in SeqIO.parse(io.StringIO(download(url)), format='fasta'):
            assert isinstance(seq, SeqRecord)
            seq.name = name  # re.sub(r"tr\|\w+\|", "", seq.name)
            yield seq


def sim(x, y):
    blosum = read_blosum()
    return blosum.loc[x, y]


def crossmap(ax1: plt.Axes, seq0: SeqRecord, seq1: SeqRecord):
    # # Make a landscape image
    # (seq0, seq1) = sorted((seq0, seq1), key=(lambda seq: len(seq.seq)))

    M = np.asarray([[sim(x, y) for y in seq1.seq] for x in seq0.seq])
    assert (M.shape == (len(seq0.seq), len(seq1.seq)))

    origin = ['upper', 'lower'][0]

    ax1.imshow(M, origin=origin)

    assert (abs(np.diff(ax1.get_ylim())) == M.shape[0])
    ax1.set_ylabel(seq0.name)

    assert (abs(np.diff(ax1.get_xlim())) == M.shape[1])
    ax1.set_xlabel(seq1.name)
    if (origin == 'upper'):
        ax1.xaxis.set_label_position('top')

    ax1.set_aspect('equal')
    ax1.set_xticks([])
    ax1.set_yticks([])


def main():
    seqs = sorted(get_seq(), key=(lambda seq: seq.id))

    print("Loaded sequences:", *[seq.name for seq in seqs], sep='\n')

    for (X, Y) in product(seqs, seqs):
        fig: plt.Figure
        ax1: plt.Axes
        (fig, ax1) = plt.subplots(figsize=(8, 8))
        crossmap(ax1, X, Y)

        filename = PARAM['figdest'] / F"{X.name.replace(' ', '_')}__vs__{Y.name.replace(' ', '_')}.png"
        filename.parent.mkdir(exist_ok=True, parents=True)
        fig.savefig(filename, **PARAM['savefig'])
        plt.close(fig)


if __name__ == '__main__':
    main()
