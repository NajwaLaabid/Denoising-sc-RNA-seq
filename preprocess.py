import pandas as pd
import anndata as ad
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt


'''
  Takes path to the given assignment data and returns a preprocessed AnnData object.
'''
def get_genes_data(path):
  df = pd.read_csv(path)
  df = df.set_index('Cell')
  # drop all missing genes
  df = df.dropna(axis=1)
  genes = ad.AnnData(X=df)

  # get cell types
  cell_types = map(lambda x: x.split('.')[0], genes.obs.index)
  set(cell_types) # we have 3 types

  arr = genes.obs.index.str.split('.').to_numpy()
  genes.obs['cells'] = np.vstack(arr.ravel())[:,1]
  genes.obs['Group'] = np.vstack(arr.ravel())[:,0]

  return genes
  
  
'''
  Performs cell and log normalization and computes pca for given anndata object.
'''
def normalize(data):
    data_norm = data.copy()
    sc.pp.normalize_per_cell(data_norm)
    sc.pp.log1p(data_norm)
    sc.pp.pca(data_norm)
    
    return data_norm
    
'''
  Generates series of plots with embedded cells.
'''
def plot(adatas, adata_labels, save=False):
    fig, axs = plt.subplots(1, len(adatas), figsize=(14,4))

    for i, (lbl, ad, ax) in enumerate(zip(adata_labels, adatas, axs)):
        sc.pl.pca_scatter(ad, color='Group', size=20, title=lbl, ax=ax, show=False, legend_loc='none')
        if i!=0: 
            ax.set_xlabel('')
            ax.set_ylabel('')
            
    plt.tight_layout()
    if save: plt.savefig('two-group-pca.pdf')