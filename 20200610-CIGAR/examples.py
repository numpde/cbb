# RA, 2020-06-14

"""
Examples of alignments from [1, p.82]

[1] Clote & Backofen, Computational Molecular Biology: An Introduction, 2000, Wiley
"""

from cigar import align, visualize

from more_itertools import first

temp = "ACACGGTCCTAATAATGGCC"
read = "ACGTACGT"
(x, y, z) = visualize(temp, read, first(align(temp, read)))
print(x, y, z, sep='\n')

temp = "CAGGAAGATCTTAGTTC"
read = "ACGTACGT"
(x, y, z) = visualize(temp, read, first(align(temp, read)))
print(x, y, z, sep='\n')
