# RA, 2020-09-12

# A look at the GSE60880:
#   Human Lung Fibroblasts treated with TGFbeta, IL1, EGF and small molecule inhibitors of TGFBR1 and p38.
#   HLF cells were incubated with each agent for 0.5 hours, 1 hour, 2 hours or 8 hours in standard conditions
#   before total RNA was prepared and hybridised to Affymetrix microarrays.
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60880
#
# Bradley G, Barrett SJ.
# CausalR: extracting mechanistic sense from genome scale data.
# Bioinformatics 2017 Nov 15;33(22):3670-3672

# Note:
# The values from the GSE "Series" are ~ those of 'expresso'
# (Goes from raw probe intensities to expression values):
# https://www.rdocumentation.org/packages/affy/versions/1.50.0/topics/expresso
# https://gist.github.com/numpde/772cd596fb5fe6036f7e29736bd1cf15

# Note:
# Potentially useful slides
# https://bioinformatics.mdanderson.org/MicroarrayCourse/Lectures/

import re, gzip
import pandas as pd
from tcga.utils import download

# Default download directory
download = download.to(rel_path="download")

url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60880/matrix/GSE60880_series_matrix.txt.gz"
with download(url).now.open(mode='rb') as gz:
    gz.seek(0)
    with gzip.open(gz, mode='r') as fd:
        sample_title = [
            re.findall(r'"([.\w]+)"', line)
            for line in fd.read().decode().splitlines()
            if line.lower().startswith("!sample_title")
        ].pop()

    gz.seek(0)  # !
    df_expr = pd.read_csv(gz, compression="gzip", comment='!', sep='\t', index_col='ID_REF').sort_index()

    assert (len(sample_title) == len(df_expr.columns))
    df_expr.columns = sample_title

# Affymetrix platform info (affyID -> gene names, etc.)
url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?mode=raw&is_datatable=true&acc=GPL571&id=17391&db=GeoDb_blob144"
with download(url).now.open(mode='r') as fd:
    df_genes = pd.read_csv(fd, sep='\t', comment='#', index_col='ID').sort_index()

assert (len(df_expr) == len(df_genes))
assert all(df_expr.index == df_genes.index)

# Index the expression array by gene symbol (instead of affyID)
df_expr = df_expr.set_index(df_genes['Gene Symbol'])

# Deal with duplicate gene symbols
# Cf., e.g.:
# Interpretation of multiple probe sets mapping to the same gene in Affymetrix GeneChips
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1784106/pdf/1471-2105-8-13.pdf
df_expr = df_expr.groupby(df_expr.index).median()
# print(list(df_expr.columns))


df_design = pd.DataFrame(
    data=[
        [c] + list(re.fullmatch(r"^(.*)_([.0-9]+)[h]*_([0-9]+)$", c).groups())
        for c in df_expr.columns
    ],
    columns=["sample", "agent", "hour", "n"],
).set_index("sample").astype({'hour': float, 'n': int})
# print(df_design)

# Sort by hour, so as to show temporal progression as a line in the plot
df_design = df_design.sort_values(by="hour")

# print(*sorted(list(df_design.index)), sep='\n')
# exit()

# print(set(df_design.agent))
print(set(df_design.hour))
# print(set(df_design.n))

import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
n_tsne_comp = 2
y = TSNE(n_components=n_tsne_comp).fit_transform(df_expr.T)
assert (y.shape == (len(df_expr.columns), n_tsne_comp))

print(set(df_design.agent))

# Note:
# - EGF (https://en.wikipedia.org/wiki/Epidermal_growth_factor)
# - TGFBeta (https://en.wikipedia.org/wiki/TGF_beta_1):
#   TGF-β is a multifunctional set of peptides that controls proliferation, differentiation, and other functions in many cell types.
# - IL1 (https://en.wikipedia.org/wiki/Interleukin-1_family)
# - GW849825X ??
#   From another report [https://journals.physiology.org/doi/prev/20131119-aop/pdf/10.1152/ajpendo.00408.2013]:
#   "The ALK4/5 inhibitor compound, GW849825X was synthesized as described in a published patent.."
#   From [https://en.wikipedia.org/wiki/ACVR1B]: ALK-4 acts as a transducer of activin or activin-like ligands (e.g., inhibin) signals.
# - GSK258899A ??
# - Cntl
#   https://www.uniprot.org/uniprot/Q9HUX4 (L-histidine 2-aminobutanoyltransferase, Pseudomonas aeruginosa)
#   https://www.uniprot.org/uniprot/A0A0D7MUC4 (Nicotianamine synthase protein, Pseudomonas aeruginosa)

agents = sorted(set(df_design.agent))
agent_color = {a: "C{}".format(n) for (n, a) in enumerate(agents)}
print(agent_color)

for ((n, agent), g) in df_design.groupby(by=['n', 'agent']):
    if agent.startswith("GSK") or agent.startswith("GW8"):
        # Unclear format of the agent label
        continue
    print(agent)
    ii = (n == df_design.n) & (agent == df_design.agent)
    (xx, yy) = (y[ii, 0], y[ii, 1])
    plt.scatter(xx, yy, s=(10 * df_design.hour[ii]), c=agent_color[agent])
    plt.plot(xx, yy, c=agent_color[agent], lw=0.5)

plt.show()
