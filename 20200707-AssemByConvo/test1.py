# RA, 2020-07-07

"""
Genome assembly via short template convolution.
Initial experiments.
"""

import io
import re
import gzip
import numpy as np
import pandas as pd

from itertools import chain
from collections import Counter

from Bio import SeqIO, SeqRecord
from tcga.utils import download
from tcga.utils import First, join
from tcga.strings import triplets

from plox import Plox

download = download.to(rel_path="../20200608-Downloads/cache")


class PARAM:
    class DATA:
        I = download(
            "ftp://ftp.ensembl.org/pub/release-100/fasta/caenorhabditis_elegans/"
            "dna/Caenorhabditis_elegans.WBcel235.dna.chromosome.I.fa.gz"
        ).now

        cDNA = download(
            "ftp://ftp.ensembl.org/pub/release-100/fasta/caenorhabditis_elegans/"
            "cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz"
        ).now


with gzip.open(io.BytesIO(PARAM.DATA.I.bytes)) as fd:
    seq = str(SeqIO.read(io.TextIOWrapper(fd), format='fasta').seq)
    assert len(seq) == 15_072_434

# with gzip.open(io.BytesIO(PARAM.DATA.cDNA.bytes)) as fd:
#     all_cdna = SeqIO.parse(io.TextIOWrapper(fd), format='fasta')
#     rec: SeqRecord
#
#     # 'F58G1.1.1 cdna chromosome:WBcel235:II:12919551:12923377:1 gene:WBGene00010263.1 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:wago-4 description:Piwi-like protein  [Source:UniProtKB/TrEMBL;Acc:O62275]'
#     p = re.compile(
#         r"(.*)chromosome:\w+:(?P<c>\w+):(?P<a>[0-9]+):(?P<b>[0-9]+):(?P<s>1|-1) gene:([.\w]+) gene_biotype:(?P<t>\w+) (.*)")
#
#     genes = pd.DataFrame([
#         p.fullmatch(rec.description).groupdict()
#         for rec in all_cdna
#     ]).astype({'a': int, 'b': int, 's': int})
#
#     # HACK
#     genes = genes[genes.s == 1]

# DEBUG
seq = seq[0:1_000_000]


# seq = seq[0:10000]


def all_kmers(string: str, k: int) -> dict:
    assert (k >= 1)
    from collections import defaultdict
    kmers = defaultdict(list)
    for i in range(0, len(string) - k + 1):
        kmer = string[i:(i + k)]
        # assert (len(kmer) == k)
        kmers[kmer].append(i)
    return dict(kmers)


k = 8
kmers = pd.Series(all_kmers(seq, k))

# Number of occurrences for each kmer
counts = kmers.apply(len).sort_values(ascending=False)
print(F"Number of {k}-mers: {len(counts)}")
print(F"Most common {k}-mers:", counts.nlargest(n=20), sep='\n')

# # kmer counts
# with Plox() as px:
#     px.a.loglog(counts.values)
#     px.show()

# # kmer distribution over the chromosome
# with Plox() as px:
#     # Number of most common kmers
#     n = 400
#
#     for (r, locs) in enumerate(kmers[counts.nlargest(n=n).index], start=1):
#         px.a.scatter(locs, r * np.ones(len(locs)), marker='o', c='k', edgecolors='none', s=0.1, alpha=1)
#
#     # # px.a.scatter(genes[genes.c == 'I'].a.values, np.zeros((genes.c == 'I').sum()), marker='o', c='g', edgecolors='none', s=0.1)
#     #
#     # for (i, gene) in genes[genes.c == 'I'].iterrows():
#     #     px.a.plot([gene.a, gene.a], [1, n], lw=0.1, c='g')
#     #     px.a.plot([gene.b, gene.b], [1, n], lw=0.1, c='r')
#
#     px.show()


# # Rarify kmers
# kmers = kmers[counts < counts.quantile(q=0.9)]

# Simulate a healthy read from the reference
read_len = 150
read = seq[100_000:][:read_len]

# Introduce an INSERT in the read
read = read[0:(read_len // 3)] + "X" + read[(read_len // 3):]

# # read location via its rarest kmers
# with Plox() as px:
#     from scipy.stats import gaussian_kde as kde
#
#     # All kmers of the read
#     read_kmers = pd.Series(all_kmers(read, k))
#
#     # The rarest ones among those
#     kmer_locs = sorted(
#         kmers[read_kmers.index],
#         key=len
#     )[0:5]
#
#     # Indicate the location of those in the genome
#     for locs in kmer_locs:
#         xx = np.arange(1, len(seq) + 1, step=100, dtype=float)
#         yy = kde(locs, bw_method=(1 / 330))(xx)
#         px.a.plot(xx, yy, '-', lw=0.2)
#     px.show()


# Brute-force convolution
with Plox() as px:
    convo = [
        (
                sum((a == b) for (a, b) in zip(seq[i:], read))
                - (1 / 4) * len(read)
        )
        for i in range(0, len(seq) - len(read) + 1)
    ]

    px.a.plot(convo, lw=0.1)
    px.show()
