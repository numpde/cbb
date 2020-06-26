# RA, 2020-06-22

"""
RNA secondary structure prediction, first algorithm of [2].
Refs: [1, p.208], [2-4].

[1] Clote & Backofen, Computational Molecular Biology: An Introduction, 2000, Wiley
[2] Nussinov et al., Algorithms for loop matching, SIAM, 1978
[3] Nussinov and Jacobson, Fast algorithm for [..] RNA, PNAS, 1980
[4] Iwakiri and Asai, RNA structure prediction, Enc of B&CB, 2018
"""

import os

import numpy as np
import pandas as pd

from itertools import product
from collections import defaultdict
from tcga.complements import dna_to_dna


def is_basepair(x, y):
    return dna_to_dna(x) == y


# class MaxAcc:
#     def __init__(self, initial=-np.inf, ties=True):
#         # self.val  # Undefined
#         # self.new  # Undefined
#         self.max = initial
#         self._ties = ties
#         self._initial = initial
#
#     def __call__(self, x):
#         self.val = x
#         self.new = (x > self.max) or ((x >= self.max) and self._ties)
#         self.max = (x if self.new else self.max)
#         return self
#
#     def __bool__(self):
#         raise RuntimeError("Please use .val, .new or .max")


# class Twain:
#     def __init__(self, s: str, minlen=1):
#         self.s = s
#         self.minlen = minlen
#
#     def __iter__(self):
#         self._i = iter(range(self.minlen, len(self.s) - self.minlen))
#         return self
#
#     def __next__(self):
#         n = next(self._i)
#         return (self.s[:n], self.s[n:])


def max_bp(s, memory={}):
    """
    Maximize number of base pairs. Refs:
    [2], [3, First algorithm], [4, Sec 'DP for MFE' on p.5]
    """

    if (len(s) <= 1):
        return {(0, "." * len(s))}

    if s not in memory:
        # Running maximum score
        m = 0
        # Candidates by score
        candidates = defaultdict(list)

        if is_basepair(s[0], s[-1]):
            for (c, p) in max_bp(s[1:-1]):
                if (c + 1 >= m):
                    m = c + 1
                    candidates[m].append("(" + p + ")")

        for j in range(1, len(s) - 1):
            for ((c1, p1), (c2, p2)) in product(max_bp(s[:j]), max_bp(s[j:])):
                if (c1 + c2 >= m):
                    m = c1 + c2
                    candidates[m].append(p1 + p2)

        # Use a set to remove duplicates
        memory[s] = {(m, p) for p in candidates[m]}

    return memory[s]


def main():
    # Example from [3, Table 1 / Fig 4]
    S = "CGGGCGGCCGGCCCCGGGCCGCGGC"
    # S = "CGGGCGGCCGGCCCCGGGCCGCG"
    # S = "CGGGCGGCCGGCCCCGGGCCG"
    # S = "CGGGCGGCCGGCCCCGGGC"
    # S = "CGGGCGGCCGGCCCCG"
    # S = "CGGGCGGCCGGCC"
    # S = "CGGGCGGC"

    r = pd.DataFrame(data=max_bp(S), columns=["score", "pairing"])
    print("Template string:", S)
    print(r)


if __name__ == '__main__':
    os.nice(19)
    main()
