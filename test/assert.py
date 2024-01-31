import anndata as ann
assert ann.read_h5ad("adata_bulk_star.h5ad").X.sum() == 73384
print('test passed')
