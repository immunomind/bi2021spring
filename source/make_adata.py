import argparse
import scanpy as sc
import warnings
import os
import pandas as pd
import numpy as np
import time


warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser()
parser.add_argument("--datadir")
parser.add_argument("--data")
parser.add_argument("--sample")
parser.add_argument("--batch")
parser.add_argument("--outadata")
parser.add_argument("--celltype")
args = parser.parse_args()
print(args.celltype)

adata1 = pd.read_csv(os.path.join(args.datadir, args.data), sep="\t", header=0, index_col=0)
adata = sc.AnnData(np.transpose(adata1))
print(adata)
print(adata.obs_names[1:3].values)
print(adata.var_names[1:3].values)
sample_adata = pd.read_csv(os.path.join(args.datadir, args.sample), header=0, index_col=0, sep="\t")
print(sample_adata.values.shape)
print(sample_adata.keys().values)
adata.obs[args.celltype] = sample_adata.loc[adata.obs_names, [args.celltype]]
adata.obs[args.batch] = sample_adata.loc[adata.obs_names, [args.batch]]
curr_time = time.time()
print(f"Saving {args.outadata}.h5ad...")
adata.write_h5ad(os.path.join(args.datadir, f"{args.outadata}.h5ad"))
print(f"Taken {time.time() - curr_time}")
