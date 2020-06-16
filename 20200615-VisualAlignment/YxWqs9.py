# RA, 2020-06-16

"""
Replicate the target–template sequence similarity in Swiss-Model.

From https://swissmodel.expasy.org/docs/help :

Target–template sequence similarity is calculated
from a normalised BLOSUM62 (Henikoff et al.) substitution matrix
(i.e. the largest and smallest values in the BLOSUM62 are 1 and 0, respectively).
The sequence similarity of the alignment is calculated
as the sum of the substitution scores divided by the number of aligned residue pairs.
Gaps are not taken into account.

Test project:
https://swissmodel.expasy.org/interactive/YxWqs9/templates/
"""

import io
import pandas as pd
import numpy as np

# Score reported by Swiss-Model
reported_score = "0.36"

join = (lambda ss: lambda j='': j.join(s.strip() for s in ss.split('\n')))

A = join("""
    QKHVVGSRSESIGGVRSAEVNTSRKRKLISDSAAEVSATVVLPVNSISTPMKWKSPRRCAVSIPKTSDEEIKEDS
    NEKLENPVISVCLEVKSKWNPKDDEQMKAVKEALHVSKAPSTVVCREDEQRRVFEFVKGCMEQKKAGSLYICGCP
    GTGKSLSMEKV--------------------------RLQAEEWAK--QAGLHCPETVSVNCTSLTKSTDIFSKI
    LGNYESGKKANGSFSPLQQLQRLFSQKQQQSRSKMMLIIADEMDYLI------TRDRGVLHELFMLTTLPLSRCI
    LIGTVFCVINVHFLKSVSYGQTSFKFKVRICPPGVANAIDLADRFLPKLK-SLNCKPLVVTFRAYSKDQILRILQ
    ERLVALPFVAFQSNALEICARKVSAASGDMRKALCVCRSALEILEIE-------VRGSID--QEPKGPV------
    --PECQ-VVKMDHMIAALSKTFKSPIVDT-IQSLPQHQQIIVCSAAKAFRGSKKDRTIAELNKLYLE-ICKSSMI
    TPAGITEFSNMCTVLNDQGILKLSLA
""")()

B = join("""
    ------------------------------------------------------------------------ESS
    PEKLQFGSQSIFLRTKA--------LLQKSSELVNLNSSDGALPARTAEYEQVMNFLAKAISEHRSDSLYITGPP
    GTGKTAQLDMIIRQKFQSLPLSLSTPRSKDVLRHTNPNLQNLSWFELPDGRLESVAVTSINCISLGEPSSIFQKI
    ---FDSFQDLNGPTLQIKNMQHLQKFLEPYHKKTTFVVVLDEMDRLLHANTSETQSVRTILELFLLAKLPTVSFV
    LIG-------------------------------MANSLDMKDRFLSRLNLDRGLLPQTIVFQPYTAEQMYEIVI
    QKMSSLPTIIFQPMAIKFAAKKCAGNTGDLRKLFDVLRGSIEIYELEKRFLLSPTRGSLNSAQVPLTPTTSPVKK
    SYPEPQGKIGLNYIAKVFSKFVNNNSTRTRIAKLNIQQKLILCTIIQSLKLN-SDATIDESFDHYIKAITKTDTL
    APLQRNEFLEICTILETCGLVSI---
""")()


def get_blosum() -> pd.DataFrame:
    # Source (2020-06-15):
    # https://raw.githubusercontent.com/wrpearson/fasta36/7c0dba1dfe5fc92d937f2bd5f9c90b8bfdb14743/data/blosum62.mat
    blosum62 = join("""
        0 A R N D C Q E G H I L K M F P S T W Y V B Z X
        A 4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0 -2 -1 0
        R -1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3 -1 0 -1
        N -2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3 3 0 -1
        D -2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3 4 1 -1
        C 0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2
        Q -1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2 0 3 -1
        E -1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2 1 4 -1
        G 0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3 -1 -2 -1
        H -2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3 0 0 -1
        I -1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3 -3 -3 -1
        L -1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1 -4 -3 -1
        K -1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2 0 1 -1
        M -1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1 -3 -1 -1
        F -2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1 -3 -3 -1
        P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2 -2 -1 -2
        S 1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2 0 0 0
        T 0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0 -1 -1 0
        W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3 -4 -3 -2
        Y -2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1 -3 -2 -1
        V 0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4 -3 -2 -1
        B -2 -1 3 4 -3 0 1 -1 0 -3 -4 0 -3 -3 -2 0 -1 -4 -3 -3 4 1 -1
        Z -1 0 0 1 -3 3 4 -2 0 -3 -3 1 -1 -3 -1 0 -1 -3 -2 -2 1 4 -1
        X 0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2 0 0 -2 -1 -1 -1 -1 -1
    """)('\n')

    blosum = pd.read_csv(io.StringIO(blosum62), sep=' ', index_col=0).astype(int)
    assert np.array_equal(blosum, blosum.T)
    return blosum


# Score matrix
S = get_blosum()

# Normalize
(min, max) = (S.min().min(), S.max().max())
S = (S - min) / (max - min)

assert (len(A) == len(B))

score = pd.Series(
    S[a][b]
    for (a, b) in zip(A, B)
    if (a in S) and (b in S)
).mean()

print("Reported score:", reported_score)
print("Computed score:", score)
# Reported score: 0.36
# Computed score: 0.3664902998236331
