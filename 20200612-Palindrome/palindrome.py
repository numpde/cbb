# RA, 2020-06-12

"""
Exercise [1, p.22]:
Find genomic palindromes in the genome of "M. jannaschii".

The first archaeon to have its complete genome sequenced.

[1] Clote & Backofen, Computational Molecular Biology: An Introduction, 2000, Wiley
[2] https://en.wikipedia.org/wiki/Methanocaldococcus_jannaschii
[3] https://science.sciencemag.org/content/273/5278/1058.long
[4] https://www.genome.jp/kegg-bin/show_organism?org=mja
[5] ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/091/665/GCA_000091665.1_ASM9166v1
[6] https://www.ncbi.nlm.nih.gov/nuccore/L77117.1
[7] https://en.wikipedia.org/wiki/FASTA_format#Sequence_representation
"""

import os
import wget
import gzip
import hashlib
import pandas as pd

from pathlib import Path
from datetime import datetime, timezone
from inclusive import range
from collections import defaultdict
from more_itertools import first

PARAM = {
    'fasta': "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/091/665/GCA_000091665.1_ASM9166v1/GCA_000091665.1_ASM9166v1_genomic.fna.gz",
    'readme': "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/091/665/GCA_000091665.1_ASM9166v1/README.txt",
    'md5': "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/091/665/GCA_000091665.1_ASM9166v1/md5checksums.txt",

    'local': Path(__file__).parent / "genome/UV/",
}

fasta_complement = {
    # Nucleotides
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',

    # pYrimidines -- puRine
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


def download(url: str):
    out = PARAM['local'] / Path(url).name
    out.parent.mkdir(exist_ok=True, parents=True)

    if out.exists():
        return out

    wget.download(str(url).format(), str(out))

    with open(str(out) + "_meta.txt", 'w') as fd:
        print("Downloaded on {} from".format(datetime.now(tz=timezone.utc)), file=fd)
        print(url, file=fd)

    return out


def download_all():
    return {
        k: download(PARAM[k])
        for k in ['fasta', 'readme', 'md5']
    }


def check_md5(files):
    md5 = pd.read_csv(files['md5'], header=None, sep=' ').dropna(axis=1).set_index(0).squeeze()
    md5 = md5.apply(lambda f: Path(f).name)
    md5 = pd.Series(index=md5, data=md5.index)

    for (k, f) in files.items():
        if (f.name in md5):
            match = (md5[f.name] == hashlib.md5(f.open('rb').read()).hexdigest())
            match = ("OK" if match else "Failed")
            print(F"md5 check {f.name}: {match}")
        else:
            print(F"No md5 check for {f.name}")


    return files


def get_genome():
    files = check_md5(download_all())

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

    with gzip.open(files['fasta'], mode='rb') as fd:
        # Get the first chunk
        (meta, genome) = first(read(fd))

    # from collections import Counter
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
    # See [7, FASTA format sequence representation]
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

    print(F"`{meta}` has length {len(genome)}")
    # L77117.1 Methanocaldococcus jannaschii DSM 2661, complete genome has length 1664970

    # kmer_histo(genome)

    counts = count_palindromes(genome)

    for (k, n) in sorted(counts.items()):
        print(F"Palindromes of length {k}: {n}")


def tests():
    assert (all_kmers("bbbc", 2) == {'bb': [0, 1], 'bc': [2]})


if __name__ == '__main__':
    os.nice(19)
    tests()
    main()
