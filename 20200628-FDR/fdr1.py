# RA, 2020-06-28

"""
Computational experiment confirming the FDR correction.
Based on [1].

[1] https://cpb-us-w2.wpmucdn.com/blog.nus.edu.sg/dist/0/3425/files/2018/10/Understanding-Benjamini-Hochberg-method-2ijolq0.pdf
"""

import numpy as np
import matplotlib.pyplot as plt

from inclusive import range
from plox import Plox
from tcga.utils import download
from pathlib import Path

# Reference [1]
slides = "https://cpb-us-w2.wpmucdn.com/blog.nus.edu.sg/dist/0/3425/files/2018/10/Understanding-Benjamini-Hochberg-method-2ijolq0.pdf"
download(slides).to(rel_path="cache/download").now


def get_obs():
    rs = np.random.RandomState(1)

    # Number of hypothesis tests
    M = 10000

    mus1 = rs.normal(size=M)
    mus2 = mus1 + (np.arange(len(mus1)) > 0.9 * len(mus1))

    # Group sizes
    s1 = 25
    s2 = 25

    obs = [
        [
            rs.normal(loc=mu1, scale=1, size=s1),
            rs.normal(loc=mu2, scale=1, size=s2),
        ]
        for (mu1, mu2) in zip(mus1, mus2)
    ]

    h0_is_true = (mus1 == mus2)

    return (obs, h0_is_true)


def pvalues(obs):
    # https://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm (http://archive.ph/bzplE)

    ts = [
        (np.mean(g1) - np.mean(g2)) / np.sqrt(np.var(g1) / len(g1) + np.var(g2) / len(g2))
        for (g1, g2) in obs
    ]

    nus = [
        # len(g1) + len(g2) - 2
        (
                (np.var(g1) / len(g1) + np.var(g2) / len(g2)) ** 2
        ) / (
                ((np.var(g1) / len(g1)) ** 2) / (len(g1) - 1)
                +
                ((np.var(g2) / len(g2)) ** 2) / (len(g2) - 1)
        )
        for (g1, g2) in obs
    ]

    from scipy.stats import t as student_t
    ps = np.asarray([
        # 2 * (0.5 - abs(0.5 - student_t.cdf(t, nu)))  # two-sided
        student_t.cdf(t, nu)  # one-sided
        for (t, nu) in zip(ts, nus)
    ])

    return ps


def main():
    (obs, h0) = get_obs()
    ps = pvalues(obs)
    ranks = np.asarray(range[1, len(ps)])

    (ps, h0, obs) = map(np.asarray, zip(*sorted(zip(ps, h0, obs))))

    alpha = 0.05
    print(F"There are {sum(ps <= alpha)} p-values less or equal {alpha}.")

    # The true False Discovery Rate for the corresponding p-values `ps`
    fdr = np.cumsum(h0) / ranks

    # Benjamini-Hochberg estimate: FDR -> p-value
    fdr_est = sorted(fdr[fdr > 0])
    ps_est = [max(ps[ps <= (a / len(ranks) * ranks)]) for a in fdr_est]

    fig_path = Path(__file__).parent / "figs"
    fig_path.mkdir(parents=True, exist_ok=True)

    with Plox() as px:
        px.a.plot(ps[fdr > 0], fdr[fdr > 0], label="True FDR")
        px.a.plot(ps_est, fdr_est, label="B-H FDR")
        px.a.plot([min(ps_est), max(ps_est)], [alpha, alpha], "r--", label=(F"FDR = {alpha}"))
        px.a.grid()
        px.a.set_xscale('log')
        px.a.set_yscale('log')
        px.a.set_xlabel("p-value")
        px.a.set_ylabel("FDR")
        px.a.legend()
        px.f.savefig(fig_path / "fdr.png")

    with Plox() as px:
        px.a.hist(ps[h0], bins='scott')
        px.a.set_title("H0 is true. Observed p-value.")
        px.f.savefig(fig_path / "h0_yay.png")

    with Plox() as px:
        px.a.hist(ps[~h0], bins='scott')
        px.a.set_title("H0 is false. Observed p-value.")
        px.f.savefig(fig_path / "h0_nay.png")

    with Plox() as px:
        px.a.hist(ps, bins='scott')
        px.a.set_title("Observed p-value.")
        px.f.savefig(fig_path / "h0_all.png")

    bh = alpha / len(ranks) * ranks
    p_fdr = max(ps[ps <= bh])

    with Plox() as px:
        px.a.scatter(ranks, ps, s=1, c=h0, cmap=plt.cm.get_cmap("copper"), label="Tests")
        px.a.plot([min(ranks), max(ranks)], [alpha, alpha], 'b-', label=F"p-value = {alpha}")
        px.a.plot([min(ranks), max(ranks)], [p_fdr, p_fdr], 'g-', label=F"FDR = {alpha}")
        px.a.plot(ranks, bh, 'r--', label="BH line")
        px.a.legend()
        px.a.grid()
        px.a.set_xscale('log')
        px.a.set_yscale('log')
        px.a.set_xlabel("rank")
        px.a.set_ylabel("p-value")
        px.f.savefig(fig_path / "p.png")


if __name__ == '__main__':
    main()
