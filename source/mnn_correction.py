import scanpy as sc
import warnings
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import time
import mnnpy


parser = argparse.ArgumentParser()
parser.add_argument("--adata")
parser.add_argument("--batch")
parser.add_argument("--celltype")
args = parser.parse_args()

warnings.filterwarnings("ignore")

# Visualisation
def save_image(resname):
    if not os.path.exists("./visualization"):
        os.makedirs("./visualization")
    resname = f"./visualization/{resname}.png"
    plt.savefig(resname, dpi=100)


def umapplot(adata, color_by, n_pcs=20, save_file_prefix="umap", n_neighbors=10, use_repx=False):
    if use_repx:
        sc.pp.neighbors(adata, use_rep='X')
    else:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    sc.tl.umap(adata)
    sc.pl.umap(adata, color=color_by, show=False)
    save_image(save_file_prefix)


def tsneplot(adata, color_by, n_pcs=20, perplexity=90, save_file_prefix="tsne", use_repx=False):
    if use_repx:
        sc.tl.tsne(adata, random_state=0, n_pcs=n_pcs, perplexity=perplexity, use_rep='X')
    else:
        sc.tl.tsne(adata, random_state=0, n_pcs=n_pcs, perplexity=perplexity)
    sc.pl.tsne(adata, color=color_by, show=False, wspace=.3)
    save_image(save_file_prefix)


# Read adata and apply some casual filters
adata = sc.read_h5ad(args.adata)
print(adata)

sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
filter_result = sc.pp.filter_genes_dispersion(adata.X, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, filter_result.gene_subset]
print(adata)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# Plot UMAP and T-SNE before correction
umapplot(adata, color_by=[args.celltype, args.batch], save_file_prefix=f"mnn_umap_{args.adata}_before_cor")

# Correction
print("Starting MNN...")
start = time.time()
adata_list = [adata[adata.obs[args.batch] == i] for i in adata.obs[args.batch].unique()]
hvgs = adata.var_names
adata_mnn = mnnpy.mnn_correct(adata_list[0], adata_list[1], var_subset=hvgs)[0]
print(f"MNN has taken {time.time() - start} seconds")

# Plot UMAP and T-SNE after correction
sc.pp.neighbors(adata_mnn, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata_mnn)
sc.pl.umap(adata_mnn, color=[args.celltype, args.batch], show=False)
resname = f"./visualization/mnn_umap_{args.adata}_after_cor.png"
plt.savefig(resname, dpi=100)

# Save corrected adata
if not os.path.exists(f"./{args.adata[:6]}"):
        os.makedirs(f"./{args.adata[:6]}")
adata_mnn.write_h5ad(os.path.join(f"./{args.adata[:6]}/mnn_{args.adata}"))
