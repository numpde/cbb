# RA, 2020-06-10

"""
Sequence alignment with CIGAR-like annotation
"""

from typing import Tuple


def align(temp: str, read: str) -> Tuple[str, int]:
    from itertools import chain, groupby
    from collections import defaultdict

    A = defaultdict(None)  # Action taken
    B = defaultdict(None)  # Backtracking
    C = defaultdict(int)  # Cost

    (visited, seenext) = ({None}, {(0, 0)})

    def propose(p):
        # i = position in the reference string
        # j = position in the read
        (i, j) = p
        (I, J) = (len(temp), len(read))
        if (i < I):
            yield ('D', (i + 1, j + 0), C[p] + 1)
        if (j < J):
            yield ('I', (i + 0, j + 1), C[p] + 1)
        if (i < I) and (j < J):
            if (temp[i] == read[j]):
                yield ('=', (i + 1, j + 1), C[p] - 1)
            else:
                yield ('X', (i + 1, j + 1), C[p] + 1)
        if (i < I):
            if (j == 0):
                yield ('S', (i + 1, j + 0), C[p] + 0)
            if (j >= J):
                yield ('S', (i + 1, j + 0), C[p] + 0)

    def onemove(p):
        for (a, q, c) in propose(p):
            if (q not in C) or (c < C[q]):
                A[q] = a
                B[q] = p
                C[q] = c
                yield q

    while seenext:
        seenext = set(chain.from_iterable(map(onemove, seenext)))

    def backtrack(q):
        while q in B:
            yield A[q]
            q = B[q]

    q = (len(temp), len(read))
    cost = C[q]

    cigar = "".join(
        "{}{}".format(k, len(list(g)))
        for (k, g) in groupby(reversed(list(backtrack(q))))
    )

    return (cigar, cost)


def visualize(temp, read, cigar: str):
    """
    Example:
        (x, y, z) = visualize("CTAGCTACCGTCGTGCTA", "GCAACGACGATG", "S3=2X1=2D1=1X1=2I1=2S3")
    returns
        x = "TACCTAGCTACCGTCG-TGCTAGC"
        y = "       GCAAC-GACGATG    "
        z = "SSSSSS==X==D=X==I==SSSSS"
    """
    import re
    i = j = 0
    x = y = z = ""
    for (a, n) in re.findall(r"[=XIDS][0-9]+", cigar):
        n = int(n)
        z += a * n
        if (a in '=XDS'):
            x += temp[i:(i + n)]
            i += n
        if (a in '=XI'):
            y += read[j:(j + n)]
            j += n
        if (a == 'I'):
            x += "-" * n
        if (a == 'S'):
            y += " " * n
        if (a == 'D'):
            y += "-" * n
    return (x, y, z)


temp = "CTAGCTACCGTCGTGCTA"
read = "GCAACGACGATG"

(cigar, cost) = align(temp, read)
(x, y, z) = visualize(temp, read, cigar)
print("Reference:", x)
print("Lone read:", y)
print("Operation:", z)
print("CIGAR:    ", cigar)
print("Cost:     ", cost)
