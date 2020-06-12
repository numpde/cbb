# RA, 2020-06-12

"""
Exercise [1, p.22]:
Find genomic palindromes in the genome of "M. jannaschii".

Run 20200612-M.jannaschii/download.py to download the genome.

[1] Clote & Backofen, Computational Molecular Biology: An Introduction, 2000, Wiley
[2] https://en.wikipedia.org/wiki/FASTA_format#Sequence_representation
"""

import os
import gzip
import pandas as pd

from pathlib import Path
from inclusive import range
from collections import defaultdict
from collections import Counter
from more_itertools import first


PARAM = {
    'fasta': Path(__file__).parent.parent / "20200608-Downloads/genomes/M.jannaschii/UV/GCA_000091665.1_ASM9166v1_genomic.fna.gz",

}

fasta_complement = {
    # Nucleotides
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',

    # pYrimidine -- puRine
    'Y': 'R',
    'R': 'Y',

    # Strong interaction: C-G
    'S': 'S',

    # Nucleic acid
    'N': 'N',

    # Ketones (G, T, U) -- bases with aMino groups (A, C)
    'M': 'K',
    'K': 'M',
}


def read(fd):
    """
    Assume the file consists of chunks like this:
        > Header 1
        a
        b
    etc. Return
        [("Header 1", "ab"), ...]
    """
    chunks = [
        (header.strip(), "".join(map(str.strip, g.split('\n'))))
        for (header, g) in [
            chunk.split('\n', 1)
            for chunk in fd.read().decode().strip().split('>') if chunk
        ]
    ]
    return chunks


def get_genome():
    with gzip.open(PARAM['fasta'], mode='rb') as fd:
        # Get the first chunk
        (meta, genome) = first(read(fd))

    # print(Counter(genome))
    # Counter({'A': 573429, 'T': 568290, 'G': 264573, 'C': 258665, 'Y': 4, 'N': 3, 'R': 3, 'M': 2, 'S': 1})

    return (meta, genome)


def all_kmers(string: str, k: int) -> dict:
    assert (k >= 1)
    kmers = defaultdict(list)
    for i in range[0, len(string) - k]:
        kmer = string[i:(i + k)]
        assert (len(kmer) == k)
        kmers[kmer].append(i)
    return dict(kmers)


def kmer_histo(genome):
    counts = pd.Series({
        k: len(all_kmers(genome, k))
        for k in range[1, 24]
    }).sort_index()

    print(counts)


def backward(genome: str) -> str:
    # See [2, FASTA format sequence representation]
    return "".join(
        fasta_complement[c]
        for c in reversed(genome)
    )


def is_genomic_palindrome(s: str):
    return (s == backward(s))


def find_palindromes(genome: str):
    K = 3

    # palindromes[k] lists the positions of palindromes of length k
    palindromes = defaultdict(list)

    # Find palindromes of length at most K - 1
    for k in range[1, K - 1]:
        for (kmer, pos) in all_kmers(genome, k).items():
            if is_genomic_palindrome(kmer):
                palindromes[k].extend(pos)

    # Find palindromes of length at least K
    # by extending palindromes of length K on both ends
    # Note: consider even and odd palindromes separately
    # Includes all nested palindromes
    for L in [K, K + 1]:
        for (kmer, ii) in all_kmers(genome, L).items():
            if not is_genomic_palindrome(kmer):
                continue
            for i in ii:
                (a, b) = (i, i + len(kmer))
                while (0 <= a < b <= len(genome)) and (genome[a] == backward(genome[b - 1])):
                    assert is_genomic_palindrome(genome[a:b])
                    palindromes[b - a].append(a)
                    a -= 1
                    b += 1

    return palindromes


def count_palindromes(genome):
    palindromes = find_palindromes(genome)

    # counts[k] is the number of palindroms locations of length k
    counts = {
        k: len(ii)
        for (k, ii) in palindromes.items()
    }

    return counts


def main():
    (meta, genome) = get_genome()
    
    # https://www.genome.jp/kegg-bin/show_organism?org=mja
    print("Note: ignoring the fact that that the chromosome is CIRCULAR")

    print(F"`{meta}` has length {len(genome)}")
    # L77117.1 Methanocaldococcus jannaschii DSM 2661, complete genome has length 1664970

    print(F"Bases content (approx): ", {k: round(v / len(genome), 4) for (k, v) in Counter(genome).items()})

    # kmer_histo(genome)

    counts = count_palindromes(genome)

    for (k, n) in sorted(counts.items()):
        print(F"Palindromes of length {k}: {n}")


def tests():
    assert (all_kmers("bbbc", 2) == {'bb': [0, 1], 'bc': [2]})
    assert (backward("TACG") == "CGTA")
    assert is_genomic_palindrome("TCGA")
    print("Tests OK")


if __name__ == '__main__':
    os.nice(19)
    tests()
    main()
