# RA, 2020-06-20

"""
A take on [1, Exercise 5.4.2].

[1] Clote & Backofen, Computational Molecular Biology: An Introduction, 2000, Wiley.
"""
from itertools import product

import pandas as pd
import numpy as np

from dataclasses import dataclass
from pathlib import Path
from more_itertools import pairwise, first, last

observations = [
    tuple("AGAAAGGTCTAGTGTTTGGTGATGTATCTATAGAGGGACG"),
    tuple("GGTCCTTTCAATATCAGTTGAATATGATGTGAGTGAGTTG"),
    tuple("GGGGGGTGGGGCCTTGATAAGAAGGGCTGTCTTTTGGTAG"),
    tuple("GTACCGGTATAGAAAAGACCGGATTCGAATTAATAATAAG"),
    tuple("TATTACTTGTTCAGCGTTATAAGATTCAGGAGGAGGTGTG"),
]


def make_stochastic(a) -> np.ndarray:
    # Normalize row-wise
    a = np.diag(1 / np.sum(a, axis=1)) @ a
    assert all(np.isclose(np.sum(a, axis=1), 1))

    return a


class HMM:
    def __init__(self, seed=0):
        self.rs = np.random.RandomState(seed=seed)

        # Hidden states
        self.Q = [1, 2, 3]

        # Output states
        self.S = ['T', 'C', 'G', 'A']

        if seed is not None:
            # Transition proba of hidden states
            self.a = pd.DataFrame(
                index=self.Q, columns=self.Q,
                data=make_stochastic(self.rs.random((len(self.Q), len(self.Q))))
            )

            # Initial proba of hidden states
            self.e = pd.Series(index=self.Q, data=make_stochastic(self.rs.random((1, len(self.Q)))).squeeze())

            # Emission proba
            self.b = pd.DataFrame(
                index=self.Q, columns=self.S,
                data=make_stochastic(self.rs.random((len(self.Q), len(self.S))))
            )

        else:
            self.a = pd.DataFrame(
                index=self.Q, columns=self.Q,
                data=[[0, 1, 0], [0, 0, 1], [1, 0, 0]]
            )

            self.e = pd.Series(index=self.Q, data=[1, 0, 0])

            self.b = pd.DataFrame(
                index=self.Q, columns=self.S,
                data=[[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0.5, 0.5]],
            )

    def __iter__(self):
        self.qq = []
        self.q = self.rs.choice(a=self.Q, p=self.e)
        return self

    def __next__(self):
        self.qq.append(self.q)
        o = self.rs.choice(a=self.S, p=self.b.loc[self.q, :])
        self.q = self.rs.choice(a=self.Q, p=self.a.loc[self.q, :])
        return o

    def __mul__(self, n: int) -> pd.Series:
        """
        Generate a sequence of length n.
        """
        return pd.Series(dict(zip(range(n), self)))


def viterbi_path_nx(hmm: HMM, O: pd.Series):
    """
    Find one of the most likely sequence of hidden states
    for the sequence of observations O.
    Use the networkx library.
    """

    import networkx as nx
    g = nx.DiGraph()

    from collections import namedtuple
    HiddenState = namedtuple('HiddenState', ['time', 'state'])

    (Alpha, Omega) = ("+", "-")

    # Negative log-likelihood
    nll = (lambda x: np.inf if (x <= 0) else -10 * np.log10(x))

    # Graph source
    for q in hmm.Q:
        # Likelihood = P[Observe the first o | Hidden state q at time 0] x P[Hidden state q at time 0]
        g.add_edge(Alpha, HiddenState(first(O.index), q), nll=nll(hmm.b[first(O)][q] * hmm.e[q]))

    # Pairwise observations
    for ((s, os), (t, ot)) in pairwise(O.items()):
        # Hidden transitions
        for (qs, qt) in product(hmm.Q, hmm.Q):
            # Likelihood = P[Observe ot | Hidden state qt at time t] x P[qt at time t | qs at time s]
            g.add_edge(HiddenState(s, qs), HiddenState(t, qt), nll=nll(hmm.b[ot][qt] * hmm.a[qt][qs]))

    # Graph sink
    for q in hmm.Q:
        g.add_edge(HiddenState(last(O.index), q), Omega, nll=0)

    # Inferred sequence of hidden states
    qq = pd.concat([
        pd.Series(index=[hs.time], data=[hs.state], dtype=object)
        for hs in nx.shortest_path(g, source=Alpha, target=Omega, weight='nll')[1:-1]
    ])

    return qq


def viterbi_path_dp(hmm: HMM, O: pd.Series):
    """
    Find one of the most likely sequence of hidden states
    for the sequence of observations O using dynamic programming.
    """

    norm = (lambda s: s / s.sum())

    d = pd.DataFrame(index=hmm.Q, columns=O.index)
    m = pd.DataFrame(index=hmm.Q, columns=O.index)
    d[first(O.index)] = hmm.b[first(O)] * hmm.e
    for ((s, __), (t, ot)) in pairwise(O.items()):
        x: pd.DataFrame
        x = hmm.a * np.outer(d[s], hmm.b[ot])
        m[t] = x.idxmax(axis=0)
        d[t] = norm(x.max(axis=0))

    # Inferred sequence of hidden states
    qq = pd.Series(index=O.index, dtype=object)

    q = d[last(d)].idxmax()
    qq[last(d)] = q
    for (s, t) in reversed(list(pairwise(O.index))):
        q = m[t][q]
        qq[s] = q

    return qq


def learn_baum_welch(hmm: HMM, O: pd.Series, niter=5, delay=0.9):
    for i in range(niter):
        # Forward variable
        a = pd.DataFrame(index=hmm.Q, columns=O.index)
        a[first(a)] = hmm.b[O[first(a)]] * hmm.e
        for (s, t) in pairwise(a):
            a[t] = hmm.b[O[t]] * (a[s] @ hmm.a)

        # Baum-Welch score Pr(O|M), based on remark in [1, p.179]:
        prom = sum(a[last(a.columns)])
        assert (prom > 0)

        print(F"Model likelihood before training step #{i + 1}: {prom}")

        # Backward variable (includes the hmm.b factor)
        b = pd.DataFrame(index=hmm.Q, columns=O.index)
        b[last(b)] = hmm.b[O[last(b)]] * 1
        for (s, t) in reversed(list(pairwise(b))):
            b[s] = hmm.b[O[s]] * (hmm.a @ b[t])

        # Remark [1, p.182]:
        if not np.isclose(prom, sum(b[first(b.columns)] * hmm.e), atol=0, rtol=1e-3):
            print("ERROR:", prom, "should equal", sum(b[first(b.columns)] * hmm.e))
            exit()

        # Expected number of transitions state i -> state j [1, Claim 5.12 and p.183]:
        n = pd.Series(
            data={
                s: hmm.a * np.outer(a[s], b[t]) / prom
                for (s, t) in pairwise(O.index)
            },
        )

        # From [1, Claim 5.9 on p.181]:
        # g = a * b / prom  # Not correct with the redefinition of b
        # Use [1, Claim 5.12 and Note on p.183]:
        g = n.apply(lambda x: x.sum(axis=1)).append(a[last(a.columns)] * 1 / prom, verify_integrity=True).T
        assert all(np.isclose(g.sum(axis=0), 1))
        assert all(np.isclose(n.sum().sum(axis=1), g[n.index].sum(axis=1)))
        assert all(np.isclose(g.groupby(O, axis=1).sum().sum(axis=1), g.sum(axis=1)))

        norm_rows = (lambda df: df.apply(lambda s: s / s.sum(), axis=1))
        hmm.e = delay * hmm.e + (1 - delay) * np.sum(first(n), axis=1)
        hmm.a = delay * hmm.a + (1 - delay) * norm_rows(n.sum())
        hmm.b = delay * hmm.b + (1 - delay) * norm_rows(pd.DataFrame(columns=hmm.S, data=g.groupby(O, axis=1).sum()).fillna(0))
        assert np.isclose(hmm.e.sum(), 1)
        assert all(np.isclose(hmm.a.sum(axis=1), 1))
        assert all(np.isclose(hmm.b.sum(axis=1), 1))


def example_viterbi():
    print('-' * 40)
    print("Example: Viterbi")

    hmm = HMM(seed=1)
    O = hmm * 20

    qq_nx = list(viterbi_path_nx(hmm, O))
    qq_dp = list(viterbi_path_dp(hmm, O))

    if not (qq_nx == qq_dp):
        print("Warning: non-identical results from viterbi_path_nx and viterbi_path_dp")
        print("Viterbi path (nx):", qq_nx, sep='\n')
        print("Viterbi path (dp):", qq_dp, sep='\n')

    qq = hmm.qq
    print("HMM's hidden state sequence:", qq, sep='\n')
    print("Viterbi path (nx):", qq_nx, sep='\n')
    print("Agreements:", [1 * (p == q) for (p, q) in zip(qq_nx, qq)], sep='\n')


def example_baum_welch():
    print('-' * 40)
    print("Example: Baum-Welch")

    # The true first-order time-homogeneous Markov model
    hmm0 = HMM(seed=None)
    # Estimate
    hmm1 = HMM(seed=1)

    for i in range(10):
        # Generate a sequence of observations from hmm0
        O = hmm0 * 100
        # Fit hmm1
        learn_baum_welch(hmm1, O, niter=15, delay=0.9)

    print(hmm0.e)
    print(hmm1.e)

    print(hmm0.a)
    print(hmm1.a)

    print(hmm0.b)
    print(hmm1.b)


def main():
    example_viterbi()
    example_baum_welch()


if __name__ == '__main__':
    main()
