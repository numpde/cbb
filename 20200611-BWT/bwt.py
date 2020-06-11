# RA, 2020-06-11

"""
Burrows-Wheeler transform
https://www.hpl.hp.com/techreports/Compaq-DEC/SRC-RR-124.pdf
https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform
"""

import numpy as np
import pandas as pd

from typing import Tuple
from more_itertools import circular_shifts, last


def btw_ref(S: str) -> Tuple[str, int]:
    shifts = sorted(circular_shifts(S))
    return ("".join(last(s) for s in shifts), shifts.index(tuple(S)))


def test_btw_ref():
    s = "^BANANA|"
    (b, i) = btw_ref(s)
    assert b == "BNN^AA|A"
    assert tuple(s) == sorted(circular_shifts(s))[i]


def ibwt_ref(L: str, I: int) -> str:
    T = sorted(range(len(L)), key=L.__getitem__)
    S = ""
    for __ in L:
        I = T[I]
        S += L[I]
    return S


def test_ibwt_ref():
    s = "^BANANA|"
    assert (s == ibwt_ref(*btw_ref(s)))


if __name__ == '__main__':
    test_btw_ref()
    test_ibwt_ref()
