import scanpy as sc
import warnings
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import time


parser = argparse.ArgumentParser()
parser.add_argument("--adata")
parser.add_argument("--batch")
parser.add_argument("--celltype")
args = parser.parse_args()

warnings.filterwarnings("ignore")

if not os.path.exists("./visualization"):
    os.makedirs("./visualization")
    
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
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata)
sc.pl.umap(adata, color=[args.celltype, args.batch], show=False)
resname = f"./visualization/combat_umap_{args.adata}_before_cor.png"
plt.savefig(resname, dpi=100)

# Correction
print("Starting Combat...")
start = time.time()
sc.pp.combat(adata, args.batch)
print(f"Combat has taken {time.time() - start} seconds")

# Plot UMAP and T-SNE after correction
#adata = sc.read_h5ad("adata4.h5ad")
#sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
#filter_result = sc.pp.filter_genes_dispersion(adata.X, min_mean=0.0125, max_mean=3, min_disp=0.5)
#adata = adata[:, filter_result.gene_subset]
#sc.pp.log1p(adata)
#sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata)
sc.pl.umap(adata, color=[args.celltype, args.batch], show=False)
resname = f"./visualization/combat_umap_{args.adata}_after_cor.png"
plt.savefig(resname, dpi=100)

# Save corrected adata
if not os.path.exists(f"./{args.adata[:6]}"):
        os.makedirs(f"./{args.adata[:6]}")
adata.write_h5ad(os.path.join(f"./{args.adata[:6]}/combat_{args.adata}"))
