# RA, 2020-06-12

"""
Follows
    [1, Sec 2.3.3: Simple statistical genomic analysis]
and
    [1, Sec 2.3.4: Genomic segmentation algorithm]

Run 20200612-M.jannaschii/download.py to download the genome.

[1] Clote & Backofen, Computational Molecular Biology: An Introduction, 2000, Wiley
"""

import io
import gzip
import numpy as np
import pandas as pd

from Bio import SeqRecord, SeqIO
from pathlib import Path
from inclusive import range
from collections import defaultdict
from collections import Counter
from more_itertools import first, windowed

PARAM = {
    'fasta': next(Path(__file__).parent.parent.glob("**/M.jannaschii/**/*genomic.fna.gz")),
}


def get_genome() -> SeqRecord:
    with gzip.open(PARAM['fasta'], mode='rb') as fd:
        for rec in SeqIO.parse(io.TextIOWrapper(fd), 'fasta'):
            return rec


def all_kmers(string: str, k: int) -> dict:
    assert (k >= 1)
    kmers = defaultdict(list)
    for (i, kmer) in enumerate(map(''.join, windowed(string, k))):
        kmers[kmer].append(i)
    return dict(kmers)


def rel_freq(kmers: dict):
    total = sum(map(len, kmers.values()))
    return {kmer: (len(ii) / total) for (kmer, ii) in kmers.items()}


def H(p: dict, b=2):
    """
    H(p) is the entropy of the discrete distribution `p` in base `b`,
    where p[A] is the frequency of event A.
    """
    p = np.array(list(p.values()))
    p = p / p.sum()
    h = -np.dot(p[p != 0], np.log(p[p != 0])) / np.log(b)
    return h


def tests():
    assert list(map(''.join, windowed("1234", 2))) == ["12", "23", "34"]
    assert all_kmers("aaab", 2) == {"aa": [0, 1], "ab": [2]}
    return True


def sliding_entropy(sequence: str, k=2, step=1, n=(2 ** 8)):
    join = (lambda s: ''.join(c for c in s if c))
    h = {
        # assert (len(x) == n)
        # assert (x == sequence[i * step: i * step + 1024])
        i: H(rel_freq(all_kmers(subseq, k)))
        for (i, subseq) in enumerate(map(join, windowed(sequence, n=n, step=step)))
    }
    h = pd.Series(h).sort_index()
    return h


def dinucleotide_entropy(sequence):
    """
    Following [1, Sec 2.3.3, p.68]
    """

    # Di-nuclear entropy over a sliding window
    import matplotlib.pyplot as plt
    (fig, ax1) = plt.subplots()
    for s in range[5, 8]:
        ax1.plot(sliding_entropy(sequence[0:2000], k=2, n=(2 ** s)), lw=1, label=F"window: $2^{s}$")
        ax1.set_xlim(0, 1000)
        ax1.set_xlabel("Window start")
        ax1.set_ylabel("Dinuclear entropy")
    ax1.legend()
    fig.savefig("dinucleotide_entropy.png")


def segmentation1(sequence):
    """
    Genomic segmentation following [1, Sec 2.3.4, p.69]
    """

    # Purines (R): A, G
    # Pyrimidines (Y): C, T, U
    RY = {'G': "R", 'A': "R", 'C': "Y", 'T': "Y", 'U': "Y", 'N': "N", 'M': "M", 'S': "S", 'R': "R", 'Y': "Y"}
    ry = (lambda seq: ''.join(map(RY.__getitem__, seq)))

    n = 717112
    sequence = ry(sequence[0:n])

    W = dict(Counter(sequence))

    js = defaultdict(int)

    (U, V) = ({k: 0 for k in W}, dict(W))
    for (m, c) in enumerate(sequence, start=1):
        U[c] += 1
        V[c] -= 1
        if (m % 100): continue
        # Jensen-Shannon divergence
        js[m] = (m / n) * (H(W) - H(U)) + ((n - m) / n) * (H(W) - H(V))

    js = pd.Series(js).sort_index()

    import matplotlib.pyplot as plt
    (fig, ax1) = plt.subplots()
    ax1.plot(js, lw=1)
    ax1.set_xlim(0, n)
    ax1.set_xlabel("Segmentation point")
    ax1.set_title(F"Jensen-Shannon divergence of the fragment {0}-{n}")
    fig.savefig("jensen-shannon.png")


def main():
    rec = get_genome()
    print(rec.description)

    # Entropies
    p1 = rel_freq(all_kmers(rec.seq, 1))
    p2 = rel_freq(all_kmers(rec.seq, 2))
    p3 = rel_freq(all_kmers(rec.seq, 3))
    print(F"H1: {H(p1)} bits")
    print(F"H2: {H(p2)} bits")
    print(F"H3: {H(p3)} bits")

    # Assuming independence
    H2_ind = H({k: p1[k[0]] * p1[k[1]] for k in p2})
    print(F"H2_ind: {H2_ind} = H2 + {H2_ind - H(p2)} bits")

    segmentation1(rec.seq)

    dinucleotide_entropy(rec.seq)


if __name__ == '__main__':
    assert tests()
    main()
