import scanpy as sc
import matplotlib.pyplot as plt
import pickle

def plot(adatas, adata_labels, filename='plot', save=False):
    fig, axs = plt.subplots(1, len(adatas), figsize=(14,4))

    for i, (lbl, ad, ax) in enumerate(zip(adata_labels, adatas, axs)):
        sc.pl.pca_scatter(ad, color='Group', size=20, title=lbl, ax=ax, show=False, legend_loc='none')
        if i!=0: 
            ax.set_xlabel('')
            ax.set_ylabel('')
            
    plt.tight_layout()
    if save: plt.savefig('img/'+filename+'.png')

'''
  Performs cell and log normalization and computes pca for given anndata object.
'''
def normalize_for_plot(data):
    data_norm = data.copy()
    sc.pp.normalize_per_cell(data_norm)
    sc.pp.log1p(data_norm)
    sc.pp.pca(data_norm)
    
    return data_norm

'''
  Takes as input the 'title of plot':'name of file' pairs in the form of a dictionary.
'''
def plot_compare_dca(save=False, filename='compare_dca', **kwargs):
    adatas = []
    adata_labels = kwargs.values()
    for d in kwargs.keys():
      with open('data/'+d+'.pickle', 'rb') as handle:
        adatas.append(normalize_for_plot(pickle.load(handle)))
    plot(adatas, adata_labels, filename=filename, save=save)