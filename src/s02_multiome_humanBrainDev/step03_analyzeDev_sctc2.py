import scanpy as sc
import numpy as np
import igraph as ig
import sctc
import matplotlib.pyplot as plt
import anndata

# Load scRNA-seq data and ensure the count matrix is a numpy array
adata = sc.read_h5ad('./rds/sctc_all.h5ad')
if not isinstance(adata.X, np.ndarray):
    adata.X = adata.X.toarray()

print(len(adata.X))

sc.pp.filter_cells(adata, min_genes=1)
sc.pp.filter_genes(adata, min_cells=1)

print(len(adata.X))

cci, gci = sctc.complexity_index(adata.X)

fo=open('./rds/sctc_all.cci','w')
for one in cci:
   fo.write(str(one)+'\n')

fo.close()



