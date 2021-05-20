import warnings


warnings.filterwarnings("ignore")

import scanpy as sc
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
import scanorama
import argparse
import bbknn
import mnnpy


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--adata")
    parser.add_argument("--batch")
    parser.add_argument("--celltype")
    parser.add_argument("--out")
    args = parser.parse_args()


def save_image(adata, resname, color_by=[args.celltype, args.batch], tool=""):
    if not os.path.exists("./visualization"):
        os.makedirs("./visualization")
    sc.pl.umap(adata, color=color_by)
    if "/" in resname:
        resname = resname.split("/")[-1]
    resname = os.path.join(".", "visualization", f"{tool}{resname[:6]}_corrected.png")
    plt.savefig(resname, dpi=100)


def save_adata(adata, resname, out, tool=""):
    if not os.path.exists(f"./{out}"):
        os.makedirs(f"./{out}")
    if "/" in resname:
        resname = resname.split("/")[-1]
    adata.write_h5ad(os.path.join(out, f"{tool}{resname[:6]}.h5ad"))


class ScanoramaCorrection:
    def __init__(self, adata, batch, celltype):
        self.adata = adata
        self.batch = batch
        self.celltype = celltype

    def filtration(self):
        print("Start filtration...\n")
        sc.pp.filter_cells(self.adata, min_genes=20)
        sc.pp.filter_genes(self.adata, min_cells=100)
        self.adata.obs["n_counts"] = self.adata.X.sum(axis=1)
        self.adata = self.adata[self.adata.obs["n_genes"] < 6000, :]
        sc.pp.normalize_per_cell(self.adata)
        self.adata.raw = sc.pp.log1p(self.adata, copy=True)
        filter_result = sc.pp.filter_genes_dispersion(self.adata.X, min_mean=0.01, max_mean=15, min_disp=0.5)
        self.adata = self.adata[:, filter_result.gene_subset]

    def correction(self):
        print("Start Scanorama...\n")
        start = time.time()
        adata_scanorama = self.adata.copy()
        adata_list = [adata_scanorama[adata_scanorama.obs[self.batch] == i] for i in adata_scanorama.obs[self.batch].unique()]
        corrected = scanorama.correct_scanpy(adata_list, return_dimred=True)
        corrected_merged_dge = corrected[0].concatenate(corrected[1:], join="outer")
        corrected_merged_dge.obs = adata_scanorama.obs
        self.adata = corrected_merged_dge
        print(f"Scanorama has taken {round(time.time() - start, 2)} seconds")

    def calculate_umap(self):
        print("Start UMAP calculation...\n")
        sc.pp.scale(self.adata, max_value=10)
        sc.tl.pca(self.adata)
        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=20)
        sc.tl.umap(self.adata)

    def fit(self):
        self.filtration()
        self.correction()
        self.calculate_umap()


class MnnCorrection:
    def __init__(self, adata, batch, celltype):
        self.adata = adata
        self.batch = batch
        self.celltype = celltype
        self.hvgs = 0
        
    def filtration(self):
        print("Start filtration...\n")
        sc.pp.filter_cells(self.adata, min_genes=20)
        sc.pp.filter_genes(self.adata, min_cells=100)
        self.adata.obs["n_counts"] = self.adata.X.sum(axis=1)
        self.adata = self.adata[self.adata.obs["n_genes"] < 6000, :]
        sc.pp.normalize_per_cell(self.adata)
        self.adata.raw = sc.pp.log1p(self.adata, copy=True)
        filter_result = sc.pp.filter_genes_dispersion(self.adata.X, min_mean=0.01, max_mean=15, min_disp=0.5)
        self.hvgs = self.adata.var_names[filter_result.gene_subset]
        self.adata = self.adata[:, filter_result.gene_subset]
        
    def correction(self):
        print("Start MNN...\n")
        start = time.time()
        mnn_adata = self.adata.copy()
        adata_list = [mnn_adata[mnn_adata.obs[self.batch] == i] for i in mnn_adata.obs[self.batch].unique()]
        corrected = mnnpy.mnn_correct(*adata_list, var_subset=self.hvgs)
        self.adata = corrected[0]
        print(f"MNN has taken {round(time.time() - start, 2)} seconds")
        
    def calculate_umap(self):
        print("Start UMAP calculation...\n")
        sc.pp.scale(self.adata, max_value=10)
        sc.tl.pca(self.adata)
        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=20)
        sc.tl.umap(self.adata)
        
    def fit(self):
        self.filtration()
        self.correction()
        self.calculate_umap()


class BbknnCorrection:
    def __init__(self, adata, batch, celltype):
        self.adata = adata
        self.batch = batch
        self.celltype = celltype
        
    def filtration(self):
        print("Start filtration...\n")
        sc.pp.filter_cells(self.adata, min_genes=20)
        sc.pp.filter_genes(self.adata, min_cells=100)
        self.adata.obs["n_counts"] = self.adata.X.sum(axis=1)
        self.adata = self.adata[self.adata.obs["n_genes"] < 6000, :]
        sc.pp.normalize_per_cell(self.adata)
        self.adata.raw = sc.pp.log1p(self.adata, copy=True)
        filter_result = sc.pp.filter_genes_dispersion(self.adata.X, min_mean=0.01, max_mean=15, min_disp=0.5)
        self.adata = self.adata[:, filter_result.gene_subset]
        
    def correction(self):
        print("Start BBKNN...\n")
        start = time.time()
        bbknn_adata = self.adata.copy()
        sc.pp.scale(bbknn_adata, max_value=10)
        sc.tl.pca(bbknn_adata)
        bbknn_adata = bbknn.bbknn(bbknn_adata, copy=True, neighbors_within_batch=5, trim=50, n_pcs=20, batch_key=self.batch, approx=False)
        print(f"BBKNN has taken {round(time.time() - start, 2)} seconds")
        
    def calculate_umap(self):
        print("Start UMAP calculation...\n")
        sc.pp.scale(self.adata, max_value=10)
        sc.tl.pca(self.adata)
        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=20)
        sc.tl.umap(self.adata)
        
    def fit(self):
        self.filtration()
        self.correction()
        self.calculate_umap()


class RegressOutCorrection:
    def __init__(self, adata, batch, celltype):
        self.adata = adata
        self.batch = batch
        self.celltype = celltype

    def filtration(self):
        print("Start filtration...\n")
        sc.pp.filter_cells(self.adata, min_genes=20)
        sc.pp.filter_genes(self.adata, min_cells=100)
        self.adata.obs["n_counts"] = self.adata.X.sum(axis=1)
        self.adata = self.adata[self.adata.obs["n_genes"] < 6000, :]
        sc.pp.normalize_per_cell(self.adata)
        self.adata.raw = sc.pp.log1p(self.adata, copy=True)
        filter_result = sc.pp.filter_genes_dispersion(self.adata.X, min_mean=0.01, max_mean=15, min_disp=0.5)
        self.adata = self.adata[:, filter_result.gene_subset]

    def correction(self):
        print("Start Regress out...\n")
        start = time.time()
        adata_reg = self.adata.copy()
        sc.pp.regress_out(adata_reg, [self.batch])
        self.adata = adata_reg
        print(f"Regress out has taken {round(time.time() - start, 2)} seconds")

    def calculate_umap(self):
        print("Start UMAP calculation...\n")
        sc.pp.scale(self.adata, max_value=10)
        sc.tl.pca(self.adata)
        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=20)
        sc.tl.umap(self.adata)

    def fit(self):
        self.filtration()
        self.correction()
        self.calculate_umap()


class CombatCorrection:
    def __init__(self, adata, batch, celltype):
        self.adata = adata
        self.batch = batch
        self.celltype = celltype
        
    def filtration(self):
        print("Start filtration...\n")
        sc.pp.filter_cells(self.adata, min_genes=20)
        sc.pp.filter_genes(self.adata, min_cells=100)
        self.adata.obs["n_counts"] = self.adata.X.sum(axis=1)
        self.adata = self.adata[self.adata.obs["n_genes"] < 6000, :]
        sc.pp.normalize_per_cell(self.adata)
        self.adata.raw = sc.pp.log1p(self.adata, copy=True)
        filter_result = sc.pp.filter_genes_dispersion(self.adata.X, min_mean=0.01, max_mean=15, min_disp=0.5)
        self.adata = self.adata[:, filter_result.gene_subset]
        
    def correction(self):
        print("Start Combat...\n")
        start = time.time()
        adata_combat = self.adata.copy()
        sc.pp.combat(adata_combat, key=self.batch)
        self.adata = adata_combat
        print(f"Combat has taken {round(time.time() - start, 2)} seconds")
        
    def calculate_umap(self):
        print("Start UMAP calculation...\n")
        sc.pp.scale(self.adata, max_value=10)
        sc.tl.pca(self.adata)
        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=20)
        sc.tl.umap(self.adata)
        
    def fit(self):
        self.filtration()
        self.correction()
        self.calculate_umap()


class Baseline:
    def __init__(self, adata, batch, celltype):
        self.adata = adata
        self.batch = batch
        self.celltype = celltype
        
    def filtration(self):
        print("Start filtration...\n")
        sc.pp.filter_cells(self.adata, min_genes=20)
        sc.pp.filter_genes(self.adata, min_cells=100)
        self.adata.obs["n_counts"] = self.adata.X.sum(axis=1)
        self.adata = self.adata[self.adata.obs["n_genes"] < 6000, :]
        sc.pp.normalize_per_cell(self.adata)
        self.adata.raw = sc.pp.log1p(self.adata, copy=True)
        filter_result = sc.pp.filter_genes_dispersion(self.adata.X, min_mean=0.01, max_mean=15, min_disp=0.5)
        self.adata = self.adata[:, filter_result.gene_subset]
        
    def calculate_umap(self):
        print("Start UMAP calculation...\n")
        sc.pp.scale(self.adata, max_value=10)
        sc.tl.pca(self.adata)
        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=20)
        sc.tl.umap(self.adata)
        
    def fit(self):
        self.filtration()
        self.calculate_umap()


if __name__ == "__main__":
    how = {
            "baseline_": Baseline,
            "combat_": CombatCorrection,
            "regress_out_": RegressOutCorrection,
            "scanorama_": ScanoramaCorrection,
            "mnn_": MnnCorrection,
            "bbknn_": BbknnCorrection
            }
    for tool in how.keys():

        adata = sc.read_h5ad(args.adata)
        adata = how[tool](adata, args.batch, args.celltype)
        adata.fit()
        save_image(adata.adata, args.adata, color_by=[args.celltype, args.batch], tool=tool)
        save_adata(adata.adata, args.adata, args.out, tool=tool)
    print("\nTHE END\n")
