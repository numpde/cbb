# RA, 2020-06-28

"""
Log-likelihood test statistic example
"""

import numpy as np
from pathlib import Path
from plox import Plox
from scipy.stats import norm, chi2, t as student_t


def experiment(n):
    # True proba density
    p = norm(loc=2, scale=1)
    # Generate a sample (observations)
    samples = p.rvs(n)
    # Proba inferred from observations
    p1 = norm(loc=np.mean(samples), scale=np.std(samples))
    # The log-likelihood statistic
    lam = 2 * (np.sum(np.log(p1.pdf(samples))) - np.sum(np.log(p.pdf(samples))))
    return lam


n = 10
r = 10000
T = [experiment(n=n) for _ in range(r)]

with Plox() as px:
    px.a.hist(T, bins='stone', density=True, label="Observed")
    df = 2
    px.a.plot(sorted(T), chi2(df=df).pdf(sorted(T)), '-', label=F"chi-squared (df={df})")
    px.a.set_title(F"Wills' theorem for {r} experiments of sample size {n}")
    px.a.set_xlabel("lambda")
    px.a.legend()
    px.f.savefig(Path(__file__).with_suffix(".png"))
