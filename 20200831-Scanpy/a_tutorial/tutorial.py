# RA, 2020-08-31

# https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html

# !mkdir data
# !wget http://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz -O data/pbmc3k_filtered_gene_bc_matrices.tar.gz
# !cd data; tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
# !mkdir write

import numpy as np
import pandas as pd
import scanpy as sc

# errors (0), warnings (1), info (2), hints (3)
sc.settings.verbosity = 3

sc.logging.print_header()

input_folder = "data/filtered_gene_bc_matrices/hg19"
results_file = "write/pbmc3k.h5ad"

# https://anndata.readthedocs.io/en/latest/anndata.AnnData.html
adata = sc.read_10x_mtx(input_folder, var_names="gene_symbols", cache=True)

print("Loaded data with size", adata.shape)

# ?
# sc.pl.correlation_matrix(adata, groupby="gene_symbols")

# # Genes with most total count?
# sc.pl.highest_expr_genes(adata, n_top=10)

# QC examples:
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0888-1
# https://bioinformatics.stackexchange.com/questions/3349/are-mitochondrial-genes-to-exclude-in-scrna-seq-such-as-ribosomal-genes
# https://www.biorxiv.org/content/10.1101/2020.02.20.958793v1

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Nomenclature convention:
# http://www.informatics.jax.org/mgihome/nomen/gene.shtml
#The mitochondria carry essential genes, among them many transfer RNA (tRNA) genes.
# Genes residing on the mitochondria have a prefix mt- ...
# print(adata.var_names)
adata.var['is_mt'] = adata.var_names.str.startswith('MT-')

non_qc_metrics_var = list(adata.var.columns)
non_qc_metrics_obs = list(adata.obs.columns)
sc.pp.calculate_qc_metrics(adata, qc_vars=['is_mt'], percent_top=None, log1p=False, inplace=True)
qc_metrics_var = list(set(adata.var.columns) - set(non_qc_metrics_var))
qc_metrics_obs = list(set(adata.obs.columns) - set(non_qc_metrics_obs))
print("QC metrics, var:", qc_metrics_var)
# QC metrics, var: ['mean_counts', 'total_counts', 'n_cells_by_counts', 'pct_dropout_by_counts']
print("QC metrics, obs:", qc_metrics_obs)
# QC metrics, obs: ['total_counts_is_mt', 'total_counts', 'n_genes_by_counts', 'pct_counts_is_mt']

# QC metrics:
# n_genes_by_counts = "the number of genes expressed in the count matrix"
# Note: adata is obs x var, i.e. sample x gene
assert list((adata.X != 0).sum(axis=1)) == list(adata.obs['n_genes_by_counts'])


# Display the QC metrics
qc_metrics = ['n_genes_by_counts', 'total_counts', 'pct_counts_is_mt']
# sc.pl.violin(adata, qc_metrics, jitter=0.4, multi_panel=True)

# sc.pl.scatter(adata, x='total_counts', y='pct_counts_is_mt')
# sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

sc.pp.regress_out
sc.tl.rank_genes_groups
