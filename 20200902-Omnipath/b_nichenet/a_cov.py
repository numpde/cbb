# RA, 2020-09-10

# Reproduce the SARS-CoV-2 case study from [1]

# [1]
# Türei et al.,
# Integrated intra- and intercellular signaling knowledge for multicellular omics analysis
# https://www.biorxiv.org/content/10.1101/2020.08.03.221242v2
#
# [2]
# RA, Notes on [1]:
# https://docs.google.com/document/d/1_qor03mpaB_ey6sUrSG4aGmSkW8WxPd7iplxh3e4uPc
#
# [3]
# Transcriptional response to SARS-CoV-2 infection
# https://doi.org/10.1101/2020.03.24.004655
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147507

# Note:
# Calu-3 is a human lung cancer cell line commonly used in cancer research and drug development.
# Calu-3 cells are epithelial and can act as respiratory models in preclinical applications.
# https://en.wikipedia.org/wiki/Calu-3

import io
import pandas as pd

from pathlib import Path
from tcga.utils import download

ROOT = Path(__file__).parent

PARAM = {
    # Transcriptome of response to SARS-CoV-2 infection
    'GSE CoV2': "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE147507&format=file&file=GSE147507%5FRawReadCounts%5FHuman%2Etsv%2Egz",

    # "111 expression datasets profiling the transcriptional response to a ligand"
    # https://zenodo.org/record/3260758
    'txn_response_ref': "https://zenodo.org/record/3260758/files/expression_settings.rds?download=1",

    'intercell': "https://omnipathdb.org/intercell",
}

# Setup default folder for downloads
download = download.to(abs_path=(ROOT / "UV/download"))


# Load the transcriptional response datasets
# Collected by NicheNet authors for optimization/validation
import tempfile
with tempfile.NamedTemporaryFile() as tf:
    tf.write(download(PARAM['txn_response_ref']).now.bytes)

    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    df = (robjects.r['readRDS'])(tf.name)
    assert (len(df) == 111)
    for (k, v) in df.items():
        print(k)
        # print(v)

# exit()

# # Doesn't work:
# import pyreadr
# result = pyreadr.read_r("~/tmpvr8xzlnf.rds", use_objects="GSE69639_TNF5_INS5_timeseries")
# exit()



def get_as_df(url, **csv_kwargs) -> pd.DataFrame:
    return pd.read_csv(io.StringIO(download(url).now.text), sep='\t', **csv_kwargs)


# The equivalent of
# https://github.com/saezlab/OmnipathR/blob/d72e73f509aa683678f5d84c07a96e2c28c5314c/R/import_Omnipath.R#L925
#
# From https://github.com/saezlab/pypath/blob/master/webservice.rst:
# The webservice currently recognizes 7 types of queries: interactions, enz_sub, annotations, complexes, intercell, queries and info...
# Interaction datasets: omnipath, pathwayextra, kinaseextra, ligrecextra (ligand-receptor interactions without literature reference), dorothea, tf_target, mirnatarget
def import_ligrecextra_interactions() -> pd.DataFrame:
    url = "https://omnipathdb.org/interactions/?datasets=ligrecextra&organisms=9606"
    return get_as_df(url)


# The equivalent of
# https://github.com/saezlab/OmnipathR/blob/d72e73f509aa683678f5d84c07a96e2c28c5314c/R/import_Omnipath.R#L1955
#
# From https://github.com/saezlab/pypath/blob/master/webservice.rst:
# Another query type is the intercell, providing information about the roles in inter-cellular signaling.
# E.g. if a protein is a ligand, a receptor, an extracellular matrix (ECM) component, etc.
def import_omnipath_intercell():
    url = "https://omnipathdb.org/intercell"
    return get_as_df(url)


# The equivalent of
# https://github.com/saezlab/OmnipathR/blob/d72e73f509aa683678f5d84c07a96e2c28c5314c/R/import_Omnipath.R#L985
def import_post_translational_interactions():
    # excludes ligrecextra
    url = "https://omnipathdb.org/interactions/?datasets=omnipath,pathwayextra,kinaseextra&organisms=9606"
    return get_as_df(url, dtype={'dip_url': str})


# The equivalent of
# https://github.com/saezlab/OmnipathR/blob/d72e73f509aa683678f5d84c07a96e2c28c5314c/R/import_Omnipath.R#L1054
def import_dorothea_interactions():
    url = "https://omnipathdb.org/interactions/?datasets=dorothea&organisms=9606&dorothea_levels=A,B,C"
    return get_as_df(url)


import_dorothea_interactions()
import_ligrecextra_interactions()
import_omnipath_intercell()
import_post_translational_interactions()
exit()


def load_cov2() -> pd.DataFrame:
    import zipfile
    with zipfile.ZipFile(download(url=PARAM['GSE CoV2']).now.local_file, mode='r') as zf:
        with zf.open("data") as fd:
            return pd.read_csv(fd, compression="gzip", sep='\t', index_col=0)


def normalize(df) -> pd.DataFrame:
    # TODO: normalize?
    return df

# RNA-seq data
df_cov2 = normalize(load_cov2())

# Mock replicates are marked True
calu3 = pd.Series(
    {
        c: ("mock" in c.lower())
        for c in df_cov2.columns
        if ("calu" in c.lower())
    },
    name="Mock",
).sort_index()
print(calu3)

# Calu-3 subset of samples
df_cov2 = df_cov2[calu3.index]
print(df_cov2)

# Make files for differential expression analysis with https://yanli.shinyapps.io/DEApp/
de_analysis_path = (ROOT / "de_analysis")
de_analysis_path.mkdir(exist_ok=True)
df_cov2.to_csv(de_analysis_path / "data.csv", sep='\t')
pd.DataFrame(calu3).reset_index().to_csv(
    de_analysis_path / "meta.csv",
    sep='\t', header=["Sample", "Treatment"], index=False
)

# (DESeq2 magic)

# Results
df_deseq2 = pd.read_csv(max(de_analysis_path.glob("DESeq2-DEG-res-*.txt")), sep='\t', index_col=0)
df_deseq2 = df_deseq2.sort_values('padj')
print(df_deseq2)

# From p.21 of https://www.biorxiv.org/content/10.1101/2020.08.03.221242v2.full.pdf
adj_p = 0.1
log_fold_change = 1
df_deseq2 = df_deseq2[(df_deseq2.log2FoldChange > log_fold_change) & (df_deseq2.padj < adj_p)]
df_cov2 = df_cov2.loc[df_deseq2.index]

# Get the ligands from the intercellular roles database
import io
from collections import Counter

df_ligands = pd.read_csv(io.StringIO(download(PARAM['intercell']).now.text), sep='\t')
df_ligands = df_ligands[df_ligands['category'] == "ligand"]
df_ligands = df_ligands[~df_ligands['genesymbol'].isna()]
# print(Counter(df_ligands['genesymbol']))

# Restrict DESeq analysis to ligands
# Restrict expression data to ligands
df_deseq2 = df_deseq2[df_deseq2.index.isin(df_ligands.genesymbol)]
df_cov2 = df_cov2[df_cov2.index.isin(df_ligands.genesymbol)]

print(len(df_cov2), "genes left:", list(df_cov2.index))
# 108 genes left: ['IL6', 'IL1A', 'CCL20', 'IL1B', 'CXCL10', 'ICAM1', 'CXCL3', 'CXCL11', 'IFNB1', 'CCL5', 'CXCL2', 'INHBA', 'TNF', 'IFNL1', 'CSF2', 'CX3CL1', 'STC2', 'CTGF', 'TNFSF10', 'CD274', 'IFNL2', 'IFNL3', 'CCL2', 'THBS1', 'LGALS9', 'SEMA7A', 'EDN1', 'HBEGF', 'C1QTNF1', 'CCL22', 'CFLAR', 'CSF1', 'LIF', 'CYR61', 'LTB', 'VEGFA', 'IL15RA', 'PLAUR', 'SEMA3A', 'BMP2', 'CXCL1', 'PDGFB', 'ULBP1', 'EDN2', 'VCAM1', 'CSF3', 'IL12A', 'TYMP', 'HLA-F', 'ADM2', 'EREG', 'S100P', 'EFNA1', 'LIFR', 'FGF19', 'FSTL3', 'BST2', 'TGFB2', 'JAG1', 'LAMA2', 'INHBE', 'MACC1', 'FST', 'SAA2', 'TNFSF15', 'IL32', 'LAMB3', 'IL15', 'L1CAM', 'SECTM1', 'ADM', 'FLT3LG', 'FGF2', 'CXCL5', 'GDF15', 'COL16A1', 'PDGFA', 'TNFSF13B', 'CXCL9', 'IL23A', 'ANGPTL4', 'IFNL4', 'PSG4', 'BDNF', 'CEACAM5', 'TNFRSF1B', 'CP', 'FLRT3', 'FAS', 'TNFSF14', 'SAA4', 'ICAM4', 'CCL28', 'VEGFC', 'PDCD1LG2', 'TNFRSF25', 'COLEC10', 'KIT', 'RLN2', 'ICAM5', 'GUCA1B', 'SEMA3D', 'DLL1', 'EBI3', 'SCUBE1', 'COL11A2', 'CYP11A1', 'IL11']
# Note: In [1], there are 117 genes.

# Note from [1, p.12]:
# Out of a total of 117 ligands over-expressed in SARS-CoV-2 infection according to NicheNet,
# we selected the 12 top-ranked ones for subsequent analysis (Supplementary Figure 6).
# Among them, we found
# o) various cytokines: interleukins (IL23A and IL1A), tumor necrosis factors (TNF and TNFSF13B)
# and
# o) chemokines (CXCL5, CXCL9 and CXCL10), known to be involved in the inflammatory response.
# TODO: according to NicheNet? -> see notes
# TODO: what is "top-ranked"? -> see notes

# "Prioritized ligands" from [1], see [1, Fig S6d] / [2, Fig S6d].
top_ranked_ref = [
    "CXCL9", "CXCL10", "CXCL5", "L1CAM", "ICAM4", "LAMA2", "TNFSF13B", "TNF", "NPPB", "INHBA", "IL1A", "IL23A",
]
top_ranked_ref = list(set(top_ranked_ref) & set(df_deseq2.index))
# print("Got all the top genes:", len(top_ranked_ref) == 12)

if True:
    df_cov2 = df_cov2.loc[top_ranked_ref]

if False:
    n_genes_to_keep = 12
    print(df_deseq2.to_markdown())
    top_de_genes = df_deseq2.log2FoldChange.nlargest(n=n_genes_to_keep).index
    df_cov2 = df_cov2.loc[top_de_genes]
    print(df_cov2)

# Done with this
del df_deseq2

#

# Get MSigDB gene sets
# https://github.com/numpde/genesets
genesets_url = "https://github.com/numpde/genesets/raw/53ce4ba8614d6d3ac2ca33243ea3f9f2c1f86ef5/genesets/msigdb/parsed/v7.1/genesets.json.zip"
df_genesets = pd.read_json(io.BytesIO(download(url=genesets_url).now.bytes), compression='zip')

# Focus on the "hallmark" gene sets, cf.
# p.22, top, https://www.biorxiv.org/content/10.1101/2020.08.03.221242v2.full.pdf
df_genesets = df_genesets[[c for c in df_genesets.columns if c.upper().startswith("HALLMARK")]]

# ... as a dictionary: gene set name -> gene set
genesets = df_genesets.T.symbols.to_dict()
# print("Gene sets:", genesets)

# Gene set enrichment analysis
# The gsea module produces GSEA results. (https://pypi.org/project/gseapy/)
# http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Main_Page
import gseapy

print(F"Running GSEA on {len(df_cov2)} genes")
# Rearrange columns to have SARS samples first, Mock samples second
df_cov2 = df_cov2[list(calu3.sort_values().index)]
gsea_res = gseapy.gsea(data=df_cov2, gene_sets=genesets, cls=list(~calu3), outdir=str(ROOT / "gsea"), min_size=2)

df_gsea_results = pd.DataFrame(
    data={
        (gs, info['fdr'], info['nes'], info['es'], info['pval'])
        for (gs, info) in gsea_res.results.items()
    },
    columns=["Geneset", "FDR", "Norm. score", "Score", "p-value"],
)
df_gsea_results = df_gsea_results.sort_values("Norm. score", ascending=False)
print(df_gsea_results.to_markdown())

# |    | Geneset                                    |      FDR |   Norm. score |     Score |   p-value |
# |---:|:-------------------------------------------|---------:|--------------:|----------:|----------:|
# |  4 | HALLMARK_IL6_JAK_STAT3_SIGNALING           | 0.341983 |      1.43446  |  0.875    | 0.0584046 |
# |  6 | HALLMARK_TNFA_SIGNALING_VIA_NFKB           | 0.204814 |      1.39715  |  0.785688 | 0.103362  |
# |  7 | HALLMARK_APOPTOSIS                         | 0.222143 |      1.34048  |  0.888889 | 0.122807  |
# |  2 | HALLMARK_INFLAMMATORY_RESPONSE             | 0.251476 |      1.28082  |  0.726963 | 0.175063  |
# |  1 | HALLMARK_ALLOGRAFT_REJECTION               | 0.216714 |      1.2367   |  0.75     | 0.153005  |
# |  3 | HALLMARK_INTERFERON_GAMMA_RESPONSE         | 0.357224 |      1.14856  |  0.777778 | 0.285266  |
# |  5 | HALLMARK_KRAS_SIGNALING_UP                 | 0.470115 |      1.00359  |  0.666667 | 0.501567  |
# |  0 | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 0.884328 |     -0.695128 | -0.444444 | 0.745665  |

leading_edge = gsea_res.results['HALLMARK_INFLAMMATORY_RESPONSE']['ledge_genes'].split(';')
