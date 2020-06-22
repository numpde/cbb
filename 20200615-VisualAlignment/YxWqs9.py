# RA, 2020-06-16

"""
Replicate the target-template sequence similarity in Swiss-Model.

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
from tcga.data.blosum import blosum62

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



# Score matrix
S = blosum62

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
