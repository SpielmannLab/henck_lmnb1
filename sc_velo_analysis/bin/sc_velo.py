#! /usr/bin/env python
import scvelo as scv
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd
import sys

# Get arguments
# path="/data/humangen_mouse/mmca/sox9_reg_ko"
#file_h5ad="/data/humangen_mouse/mmca/sox9_reg_ko/obj_sox9_wt_WT.h5ad"
#plot_filetype=".pdf"
#proportions_groupby="sub_trajectory"

file_h5ad=sys.argv[1]
plot_filetype=sys.argv[2]
subset_key = sys.argv[3]
subset_val = sys.argv[4]

color_by_before_split = 'condition'
split_by = 'condition'
color_by_after_split = "annotation"

print("file_h5ad="+file_h5ad)
print("plot_filetype="+plot_filetype)

# Set the base name for all plots
tmp_name=os.path.basename(file_h5ad)
plot_filename=os.path.splitext(tmp_name)[0] + \
    "_subset_" + subset_key + "_" + subset_val + plot_filetype
print(plot_filename)

# Read the h5ad file
adata = scv.read(file_h5ad, cache=True) # This is the file created by seurat_to_anndata.R code

# Subset based on the subset_key and subset_val provided
adata = adata[adata.obs[subset_key] == subset_val]

# Convert annotations and conditions to pandas categorical
adata.obs[color_by_before_split] = pd.Categorical(adata.obs[color_by_before_split]) 
adata.obs[color_by_after_split] = pd.Categorical(adata.obs[color_by_after_split]) 
adata.obs[split_by] = pd.Categorical(adata.obs[split_by]) 

sc.pl.umap(adata, color=color_by, frameon=False, save=plot_filename, legend_loc = "right margin")

## Normalize, do PCA, and caluclate veloctiy and save plot for entire dataset
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=50, n_neighbors=30)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis='umap', color=color_by_before_split, legend_loc="right margin", save=plot_filename)

# Now do the same thing for one split_by category at a time
print("Plotting one " + split_by + " at a time: ")
for condition_of_interest in adata.obs.loc[:, split_by].unique():
    print(condition_of_interest)
    adata_sub = adata[adata.obs[split_by] == condition_of_interest]
    
    scv.pp.neighbors(adata_sub, n_pcs = 50, n_neighbors = 30)
    scv.pp.moments(adata_sub, n_pcs=50, n_neighbors=30)
    scv.tl.velocity(adata_sub)
    scv.tl.velocity_graph(adata_sub)
    
    plot_filename_sub = plot_filename.strip(plot_filetype) + "_splitby_" + split_by + "_" + condition_of_interest + plot_filetype
    scv.pl.velocity_embedding_stream(adata_sub, basis='umap', color=color_by, legend_loc="right margin", save=plot_filename_sub)

