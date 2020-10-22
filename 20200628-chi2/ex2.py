# RA, 2020-10-22

"""
Log-likelihood test statistic example.
"""

import numpy as np
from pathlib import Path
from plox import Plox
from scipy.stats import norm, chi2, linregress


def experiment(n):
    # True model
    a = -1
    b = 0.0
    x = np.arange(0, n)
    s = 1
    # Generate observations
    y = a + b * x + norm(loc=0, scale=s).rvs(n)
    # Linear regression
    (b_hat, a_hat, _, _, _) = linregress(x, y)
    # The log-likelihood statistic
    lam = (s ** (-2)) * (np.sum(np.power(y - np.mean(y), 2)) - np.sum(np.power(y - (a_hat + b_hat * x), 2)))
    return lam


n = 10
r = 1000
T = np.array([experiment(n=n) for _ in range(r)])

with Plox() as px:
    T = T[(np.quantile(T, 0.01) <= T) & (T <= np.quantile(T, 0.9))]
    px.a.hist(T, bins='stone', density=True, label="Observed")
    df = 1
    px.a.plot(sorted(T), chi2(df=df).pdf(sorted(T)), '-', label=F"chi-squared (df={df})")
    px.a.set_title(F"Wilks' theorem for {r} experiments of sample size {n}")
    px.a.set_xlabel("lambda")
    px.a.legend()
    px.f.savefig(Path(__file__).with_suffix(".png"))
