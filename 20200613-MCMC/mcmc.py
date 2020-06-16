# RA, 2020-06-13

"""
MCMC for sampling a distribution
"""

import numpy as np
import pandas as pd
import typing
import matplotlib.pyplot as plt

from os import nice
from joblib import delayed, Parallel
from pathlib import Path
from inclusive import range
from contextlib import contextmanager
from progressbar import progressbar

root_path = Path(__file__).parent
savefig_param = dict(bbox_inches='tight', pad_inches=0, transparent=False)


def mcmc(scenario, seed):
    rs = np.random.RandomState(seed)

    # Maximal number of MCMC steps
    K = 2 ** 22


    # Interval domain of L
    (a, b) = (0, 1)

    # Likelihood
    def L(x):
        return 1 + np.sin(x * 5)

    assert (0 <= scenario <= 7)

    # scenario 0
    if True:
        # Number of states
        n = 11
        # States
        J = list(range[1, n])
        # Neighbors
        N = {i: [j for j in (i - 1, i + 1) if j in J] for i in J}
        # x-coordinate of states
        X = dict(zip(J, np.linspace(a, b, 2 * len(J) + 1)[1::2]))
        # MCMC step acceptance rule of moving from state i to j
        accept = (lambda i, j: rs.uniform(0, P[i]) <= P[j])

    if (scenario >= 1):
        N = {i: [(j if j in J else i) for j in (i - 1, i + 1)] for i in J}

    if (scenario >= 2):
        N = {i: [j for j in (i - 4, i - 1, i + 1, i + 4) if j in J] for i in J}
        accept = (lambda i, j: rs.uniform(0, P[i] / len(N[i])) <= P[j] / len(N[j]))

    if (scenario >= 3):
        X = dict(zip(J, sorted(rs.uniform(a, b, len(J)))))

    if (scenario >= 4):
        xx = np.linspace(a, b, 111)
        X = dict(zip(J, sorted(rs.choice(xx, size=len(J), replace=True, p=(L(xx) / sum(L(xx)))))))

    if (scenario >= 5):
        X = dict(zip(J, rs.uniform(a, b, len(J))))

    if (scenario >= 6):
        # Number of states
        n = 97
        # States
        J = list(range[1, n])
        # Neighbors
        N = {i: [j for j in (i - 4, i - 1, i + 1, i + 4) if j in J] for i in J}
        # x-coordinate of states
        X = dict(zip(J, np.linspace(a, b, 2 * len(J) + 1)[1::2]))

    if (scenario >= 7):
        # x-coordinate of states
        X = dict(zip(J, rs.uniform(a, b, len(J))))

    # Reference frequency of states
    P = {i: L(X[i]) for i in J}

    # Observed visit frequency
    Q = {i: 0 for i in J}

    # Initial state
    i = rs.choice(J)

    # Reference frequency distribution (unnormalized)
    p = pd.Series(index=[X[j] for j in J], data=[P[j] for j in J])

    @contextmanager
    def plot_dist(q, p) -> typing.Iterable[plt.Figure]:
        p = p / (p.sum() or 1)
        q = q / (q.sum() or 1)
        fig: plt.Figure
        ax1: plt.Axes
        (fig, ax1) = plt.subplots(figsize=(8, 5))
        ax1.clear()
        ax1.scatter(p.index, p, marker='o', c='C0', s=36, zorder=10, label="Reference")
        ax1.scatter(q.index, q, marker='o', c='C3', s=30, zorder=20, label="MCMC")
        xx = np.linspace(a, b, 1001)
        yy = L(xx) / sum(L(q.index))
        ax1.plot(xx, yy, '--', c='C0', lw=1)
        ax1.plot([a, b], [0, 0], lw=1, c='k', zorder=-10)
        ax1.set_xlim(a, b)
        ax1.set_xticks([])
        ax1.set_ylim(round(-0.1 * max(yy), 2), round(1.35 * max(yy), 2))
        ax1.set_yticks([])
        # ax1.axis('off')
        ax1.legend(loc='upper right')
        yield fig
        plt.close(fig)

    # Convergence history
    h = pd.Series(dtype=float, name="chi", index=pd.Series(dtype=int, name="k"))

    def chi(q: pd.Series, p: pd.Series):
        # Eqn (4) in
        # https://www.cse.huji.ac.il/~werman/Papers/ECCV2010.pdf
        p = p / p.sum()
        q = q / q.sum()
        x = np.sqrt(0.5 * (((p - q) ** 2) / (p + q)).sum())
        assert (0 <= x <= 1)
        return x

    def checkpoint():
        q = pd.Series(index=[X[j] for j in J], data=[Q[j] for j in J])
        d = chi(p, q)
        h.at[k] = d

        out_path = root_path / F"./figs_rrt/scenario={scenario:03}/seed={seed:04}"
        out_path.mkdir(exist_ok=True, parents=True)
        # print("Writing to", out_path)

        with plot_dist(q, p) as fig:
            assert isinstance(fig, plt.Figure)
            fig.savefig(out_path / F"k={k:09}.png", **savefig_param)

        h.to_csv(out_path / F"chi_hist.csv", sep='\t')

    # Metropolis–Hastings
    for k in progressbar(range[0, K]):
        if (k & (k - 1) == 0):
            checkpoint()

        # Record the visit to state i
        Q[i] += 1

        # Next state proposal
        j = rs.choice(N[i])

        # 
        if accept(i, j):
            i = j

    return h


def run_scenario(scenario, n_jobs=8):
    print(F"Scenario {scenario}")
    Parallel(n_jobs=n_jobs)(
        delayed(mcmc)(scenario, seed)
        for seed in range[1, 8]
    )


def main():
    # # Testrun
    # return run_scenario(4, n_jobs=1)

    for scenario in range[0, 7]:
        run_scenario(scenario)


if __name__ == '__main__':
    nice(19)
    main()
