import warnings

warnings.filterwarnings("ignore")
import argparse
import scanpy as sc
import os
import matplotlib.pyplot as plt
import time
import scanorama
import bbknn
import mnnpy
from scipy.stats import mannwhitneyu
from sklearn.metrics import silhouette_samples



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--adata",
                        help="""Expression matrix with n_obs x n_vars shape,
                                where: 
                                n_obs is the number of cells
                                n_vars is the number of genes.
                                Should be in h5ad format and also have "batch" and "cell type" columns.""")
    parser.add_argument("--batch",
                        help="""The name of column which reflects batch labels""")
    parser.add_argument("--celltype",
                        help="""The name of column which reflects cell type labels""")
    parser.add_argument("--out",
                        help="""Directory for saving corrected expression matrix and some useful UMAP plots""")
    parser.add_argument("--do_filter", type=bool, choices=[True, False], default=True,
                        help="""This argument tells if programm will apply some basic scRNA data analysis filters 
                             (e.g. gene with more than 3 cells etc.)""")
    parser.add_argument("--algtorun", choices=["all", "combat", "mnn", "bbknn", "scanorama", "regress_out"],
                        help="""This arguments means which algorithms program will use to correct batch effect.
                                If "all" is chosen program will run all (5) correction methods and then
                                print the name of the best one in stdout""")

    args = parser.parse_args()


def save_image(adata, resname, color_by=(args.celltype, args.batch), tool=""):
    if not os.path.exists(f"{args.out}/visualization"):
        os.makedirs(f"{args.out}/visualization")
    sc.pl.umap(adata, color=color_by)
    if "/" in resname:
        resname = resname.split("/")[-1]
    resname = os.path.join(args.out, "visualization", f"{tool}{resname[:6]}_corrected.png")
    plt.savefig(resname, dpi=100)


def save_adata(adata, resname, out, tool=""):
    if not os.path.exists(f"./{out}"):
        os.makedirs(f"./{out}")
    if "/" in resname:
        resname = resname.split("/")[-1]
    adata.write_h5ad(os.path.join(out, f"{tool}{resname[:6]}.h5ad"))


def calculate_significance(alg_adata, base_adata_silhouettes, by: str, adj_alpha: float, alternative: str):
    """
    :param alg_adata : Expression matrix, corrected with particular algorithms.
    :param base_adata_silhouettes: Silhouettes of expression matrix before correction.
    :param by: Batch or cell type labels of samples.
    :param adj_alpha: Adjusted p-value for MU-test.
    :param alternative: less or greater based if by is equal cell type or batch.
    :return: p-value of MU-test < adjusted alpha.
    """
    alg_silhouettes = silhouette_samples(alg_adata.X, labels=alg_adata.obs[by], metric="cosine")
    p_val = mannwhitneyu(base_adata_silhouettes, alg_silhouettes, alternative=alternative).pvalue

    return int(p_val < adj_alpha)


class ScanoramaCorrection:
    def __init__(self, adata, batch, celltype, do_filter=True):
        self.adata = adata
        self.batch = batch
        self.celltype = celltype
        self.do_filter = do_filter

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
        adata_list = [adata_scanorama[adata_scanorama.obs[self.batch] == i] for i in
                      adata_scanorama.obs[self.batch].unique()]
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
        if self.do_filter:
            self.filtration()
        self.correction()
        self.calculate_umap()


class MnnCorrection:
    def __init__(self, adata, batch, celltype, do_filter=True):
        self.adata = adata
        self.batch = batch
        self.celltype = celltype
        self.hvgs = 0
        self.do_filter = do_filter

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
        if self.do_filter:
            self.filtration()
        self.correction()
        self.calculate_umap()


class BbknnCorrection:
    def __init__(self, adata, batch, celltype, do_filter=True):
        self.adata = adata
        self.batch = batch
        self.celltype = celltype
        self.do_filter = do_filter

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
        bbknn_adata = bbknn.bbknn(bbknn_adata, copy=True, neighbors_within_batch=5, trim=50, n_pcs=20,
                                  batch_key=self.batch, approx=False)
        print(f"BBKNN has taken {round(time.time() - start, 2)} seconds")

    def calculate_umap(self):
        print("Start UMAP calculation...\n")
        sc.pp.scale(self.adata, max_value=10)
        sc.tl.pca(self.adata)
        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=20)
        sc.tl.umap(self.adata)

    def fit(self):
        if self.do_filter:
            self.filtration()
        self.correction()
        self.calculate_umap()


class RegressOutCorrection:
    def __init__(self, adata, batch, celltype, do_filter=True):
        self.adata = adata
        self.batch = batch
        self.celltype = celltype
        self.do_filter = do_filter

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
        if self.do_filter:
            self.filtration()
        self.correction()
        self.calculate_umap()


class CombatCorrection:
    def __init__(self, adata, batch, celltype, do_filter=True):
        self.adata = adata
        self.batch = batch
        self.celltype = celltype
        self.do_filter = do_filter

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
        if self.do_filter:
            self.filtration()
        self.correction()
        self.calculate_umap()


class Baseline:
    def __init__(self, adata, batch, celltype, do_filter=True):
        self.adata = adata
        self.batch = batch
        self.celltype = celltype
        self.do_filter = do_filter

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
        if self.do_filter:
            self.filtration()
        self.calculate_umap()


if __name__ == "__main__":
    COMPARISONS = 5 * 2  # 5 algorithms and 2 label types: batch and cell type
    ADJ_ALPHA = 0.05 / COMPARISONS

    how = {
        "combat_": CombatCorrection,
        "regress_out_": RegressOutCorrection,
        "scanorama_": ScanoramaCorrection,
        "mnn_": MnnCorrection,
        "bbknn_": BbknnCorrection
    }

    result = []

    adata = sc.read_h5ad(args.adata)
    adata = Baseline(adata, args.batch, args.celltype, do_filter=args.do_filter)
    adata.fit()
    save_image(adata.adata, args.adata, color_by=(args.celltype, args.batch), tool="baseline_")
    save_adata(adata.adata, args.adata, args.out, tool="baseline_")
    baseline_silhouettes_batch = silhouette_samples(adata.adata.X, labels=adata.adata.obs[args.batch], metric="cosine")
    baseline_silhouettes_cell = silhouette_samples(adata.adata.X, labels=adata.adata.obs[args.celltype], metric="cosine")

    if args.algtorun == "all":
        for tool in how.keys():
            adata = sc.read_h5ad(args.adata)
            adata = how[tool](adata, args.batch, args.celltype, do_filter=args.do_filter)
            adata.fit()
            save_image(adata.adata, args.adata, color_by=[args.celltype, args.batch], tool=tool[:-1])
            save_adata(adata.adata, args.adata, args.out, tool=tool[:-1])

            cell = calculate_significance(adata.adata, baseline_silhouettes_cell,
                                          by=args.celltype, adj_alpha=ADJ_ALPHA,
                                          alternative="less")
            batch = calculate_significance(adata.adata, baseline_silhouettes_batch,
                                          by=args.batch, adj_alpha=ADJ_ALPHA,
                                          alternative="greater")
            result.append((tool[:-1], cell * 1.5 + batch))

    else:
        adata = sc.read_h5ad(args.adata)
        adata = how[args.algtorun](adata, args.batch, args.celltype, do_filter=args.do_filter)
        adata.fit()
        save_image(adata.adata, args.adata, color_by=(args.celltype, args.batch), tool=f"{args.algtorun}_")
        save_adata(adata.adata, args.adata, args.out, tool=f"{args.algtorun}_")

    print("\nTHE END OF CORRECTION\n")
    result.sort(key=lambda x: x[1])
    print("The best algorithms is: ")
    for alg, res in result:
        if res >= 1.5:
            print(f"1. {alg.upper()}")
    print("We suggest you to try one of them.")

