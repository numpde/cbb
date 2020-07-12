# RA, 2020-07-05

"""
ORF detection by looking at base frequencies?

(Initial draft)
"""

import io
import gzip
import numpy as np
import pandas as pd

from collections import Counter

from Bio import SeqIO
from tcga.utils import download
from tcga.utils import First
from tcga.strings import triplets

download = download.to(rel_path="../20200608-Downloads/cache")


class PARAM:
    class DATA:
        I = download(
            "ftp://ftp.ensembl.org/pub/release-100/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.chromosome.I.fa.gz"
        ).now


with gzip.open(io.BytesIO(PARAM.DATA.I.bytes)) as fd:
    rec: SeqIO.SeqRecord
    rec = SeqIO.read(io.TextIOWrapper(fd), format='fasta')


from tcga.codons import standard_rna as rna_to_aa
from tcga.complements import dna_to_rna, dna_to_dna
f = First(dna_to_dna).then(dna_to_rna).then(triplets).each(rna_to_aa)

from plox import Plox
import matplotlib.pyplot as plt
from typing import Tuple

fig: plt.Figure
AX: Tuple[plt.Axes]
[fig, AX] = plt.subplots(3, 1)

for (offset, ax) in zip([0, 1, 2], AX):
    K = 3 * 100
    N = 100
    from tcga.strings import kplets
    df = pd.DataFrame(data=[
        dict(Counter(f(kplet)))
        for (kplet, _) in zip(kplets(K, rec.seq[offset:], no_remainder=False), range(N))
    ]).T

    totals = df.sum(axis=1).sort_values(ascending=False)
    print(F"Offset: {offset}, totals: {list(totals.items())}")

    df = df.loc[totals.index]

    ax.imshow(np.log(df))

plt.show()

