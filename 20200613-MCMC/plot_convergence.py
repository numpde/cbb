# RA, 2020-06-14

import re
import os
import pandas as pd
import matplotlib.pyplot as plt

from decimal import Decimal
from pathlib import Path
from more_itertools import first, last

root_path = Path(__file__).parent
file_pattern = "figs/scenario={scenario}/seed={seed}/chi_hist.csv"
savefig_param = dict(bbox_inches='tight', pad_inches=0, transparent=False)


def undo(pattern, string) -> dict:
    from string import Formatter
    return re.fullmatch(
        pattern.format(**{
            k: r"(?P<{}>.+)".format(k)
            for (__, k, __, __) in Formatter().parse(pattern)
        }),
        string
    ).groupdict()


all_files = pd.DataFrame(
    {'file': f, **undo(file_pattern, os.path.relpath(str(f), root_path))}
    for f in root_path.glob(file_pattern.format(scenario='*', seed='*'))
)

# Group by scenario
all_files = {
    scenario: pd.concat(axis=1, objs=[
        pd.read_csv(
            file_info.file,
            sep='\t',
            index_col=0, header=0,
            names=["k", file_info.seed],
        )
        for (__, file_info) in files.iterrows()
    ])
    for (scenario, files) in all_files.groupby('scenario')
}

(fig, ax1) = plt.subplots(figsize=(8, 5))

# Sort by final error
all_files = sorted(all_files.items(), key=(lambda __df: last(__df[1].mean(axis=1))), reverse=True)

for (n, (scenario, df)) in enumerate(all_files):
    c = F"C{n}"
    df = df[df.index > 0]
    m = df.mean(axis=1)
    s = df.std(axis=1)
    x = m.index * ((1.02) ** (n % 3))
    ax1.loglog(x, m, c=c, lw=2, ls='-', marker='o', ms=5, label=F"Scenario {Decimal(scenario)}")
    ax1.errorbar(x=x, y=m, lw=1, yerr=s, c=c)

# Use the last m
x = pd.Series(m.index).quantile([0.3, 0.9], interpolation='nearest')
ax1.loglog(x, (x ** (-1/2)), lw=1, ls='--', c='k', marker=None, label="$\mathrm{iteration}^{-1/2}$")

ax1.grid()
ax1.set_xlabel("MCMC iteration")
ax1.set_ylabel("chi error")
ax1.legend(
    # *map(reversed, ax1.get_legend_handles_labels()),
    loc='lower left', title="By final error:",
)

fig_file = root_path / Path(file_pattern).parts[0] / "convergence.png"
fig.savefig(fig_file, **savefig_param)

print(F"Saved to {fig_file}")
