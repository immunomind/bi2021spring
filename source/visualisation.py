import scanpy as sc
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def save_image(basename):
    if not os.path.exists("./visualization"):
        os.makedirs("./visualization")
    resname = f"./visualisation/{resname}.png"
    plt.savefig(resname, dpi=100)
    plt.close()


def tsneplot(adata, color_by, n_pcs=20, perplexity=90, save_file_prefix="tsne", use_repx=False):
    if use_repx:
        sc.tl.tsne(adata, random_state=0, n_pcs=n_pcs, perplexity=perplexity, use_rep='X')
    else:
        sc.tl.tsne(adata, random_state=0, n_pcs=n_pcs, perplexity=perplexity)
    sc.pl.tsne(adata, color=color_by, show=False, wspace=.3)
    save_image(save_file_prefix)


def umapplot(adata, color_by, n_pcs=20, save_file_prefix="umap", n_neighbors=10, use_repx=False):
    import scanpy as sc
    if use_repx:
        sc.pp.neighbors(adata, use_rep='X')
    else:    
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    sc.tl.umap(adata)
    sc.pl.umap(adata, color=color_by)
#    save_image(save_file_prefix)

