import scirpy as ir
import scanpy as sc
from glob import glob
import pandas as pd
import tarfile
import anndata
import warnings
import matplotlib.pyplot as plt

sc.set_figure_params(figsize=(4, 4))
sc.settings.verbosity = 4  # verbosity: errors (0), warnings (1), info (2), hints (3)

# define sample metadata. Usually read from a file.
samples = {
    "NT-lsdB": {"group": "NT-lsdB-C"},
    "LAC-alum": {"group": "LAC-alum-C"},
    "LAC-lsdB": {"group": "LAC-lsdB-2-C"}
}

# Create a list of AnnData objects (one for each sample)
adatas = []
for sample, sample_meta in samples.items():
    gex_file = "/data/CRprocessed_BCR/"+sample_meta["group"]+"_ge/filtered_feature_bc_matrix.h5"
    tcr_file = "/data/CRprocessed_BCR/"+sample+"_vdj/outs/filtered_contig_annotations.csv.gz"
    adata = sc.read_10x_h5(gex_file)
    adata_tcr = ir.io.read_10x_vdj(tcr_file)
    ir.pp.merge_with_ir(adata, adata_tcr)
    adata.obs["sample"] = sample
    adata.obs["group"] = sample_meta["group"]
    # concatenation only works with unique gene names
    adata.var_names_make_unique()
    adatas.append(adata)

# Merge anndata objects
adata = anndata.concat(adatas)

#print basic adata info
print("adata obs columns:")
print(adata.obs.columns)
print("adata shape:")
print(adata.shape)
print("adata var_names:")
print(adata.var_names)

######################################################################################################################
#data preprocessing follows the scanpy tutorial and the best practice paper by Luecken et al (https://www.embopress.org/doi/full/10.15252/msb.20188746).
#Filter genes based on number of cells or counts.
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.filter_cells(adata, min_genes=100)
#Normalize counts per cell.
sc.pp.normalize_total(adata)
#Logarithmize the data matrix.
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="cell_ranger", n_top_genes=5000)
sc.tl.pca(adata)
sc.pp.neighbors(adata)


#Compute CDR3 neighborhood graph and define clonotypes
for cutoff in [9]:
	ir.pp.ir_dist(
	    adata,
	    metric="alignment",
	    sequence="aa",
	    cutoff=cutoff,
	)
	adata.obs_names_make_unique()
	ir.tl.define_clonotype_clusters(
	    adata, sequence="aa", metric="alignment", receptor_arms="all", dual_ir="any"
	)

	for min_cells in [20]:
		ir.tl.clonotype_network(adata, min_cells=min_cells, sequence="aa", metric="alignment")
		ax = ir.pl.clonotype_network(
		    adata, 
		    color="sample",
		    size_power = 0.5
		)
		plt.savefig("results/"+str(cutoff)+"_clonotype_network_"+str(min_cells),bbox_inches='tight',dpi=100)
