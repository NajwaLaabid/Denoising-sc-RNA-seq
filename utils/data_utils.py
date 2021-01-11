from rpy2.robjects.packages import STAP
import rpy2.robjects as ro
import os
import scipy as sp
import scanpy as sc 
from sklearn.model_selection import train_test_split
import pickle
import numpy as np 
import pandas as pd

def get_simdata(name='sim'):
	data_dir = 'data/'
	if name=='sim': 
		path_raw = data_dir+name+'_raw_2.pickle'
		path_true = data_dir+name+'_true_2.pickle'
	else:
		path_raw = data_dir+name+'.pickle'
		path_true = data_dir+name+'.pickle'

	assert os.path.exists(path_raw) and os.path.exists(path_true), "The data you're requesting doesn't exist. Use generate_simdata to generate data."

	with open(path_raw, 'rb') as handle:
		raw = pickle.load(handle)

	with open(path_true, 'rb') as handle:
		true = pickle.load(handle)

	return raw, true

# function to convert rpy2 object to pandas
def convert_rpy2(r_df):
    with localconverter(ro.default_converter + pandas2ri.converter):
        df = ro.conversion.rpy2py(r_df)
    return df

def generate_simdata(name='sim', nGroups=2, nGenes=200, batchCells=2000, dropout=5):
    with open('utils/simulate.r', 'r') as f:
        string = f.read()
    simulate = STAP(string, "simulate")

    data = simulate.simulate(nGroups=nGroups, nGenes=nGenes, batchCells=batchCells, dropout=dropout)

    # data_pd = dict((k, convert_rpy2(v)) for k, v in zip(data.names, list(data)))

    # raw = sc.AnnData(data_pd['counts'].values, obs=data_pd['cellinfo'], var=data_pd['geneinfo'])
    # raw.obs_names = data_pd['cellinfo'].Cell
    # raw.var_names = data_pd['geneinfo'].Gene
    # sc.pp.filter_genes(raw, min_counts=1)
    # with open('data/'+name+'_'+str(nGroups)+'_raw.pickle', 'wb') as handle:
    #     pickle.dump(raw, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # true = sc.AnnData(data_pd['truecounts'].values, obs=data_pd['cellinfo'], var=data_pd['geneinfo'])
    # true.obs_names = data_pd['cellinfo'].Cell
    # true.var_names = data_pd['geneinfo'].Gene
    # true = true[:, raw.var_names].copy()
    # with open('data/'+name+'_'+str(nGroups)+'_true.pickle','wb') as handle:
    #     pickle.dump(true, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
    return data

def read_dataset(adata, test_split=False, copy=False):
    if copy: adata = adata.copy()

    # check if observations are unnormalized using first 10
    X_subset = adata.X[:10]
    norm_error = 'Make sure that the dataset (adata.X) contains unnormalized count data.'
    if sp.sparse.issparse(X_subset):
        assert (X_subset.astype(int) != X_subset).nnz == 0, norm_error
    else:
        assert np.all(X_subset.astype(int) == X_subset), norm_error

    if test_split:
        train_idx, test_idx = train_test_split(np.arange(adata.n_obs), test_size=0.1, random_state=42)
        spl = pd.Series(['train'] * adata.n_obs)
        spl.iloc[test_idx] = 'test'
        adata.obs['dca_split'] = spl.values
    else:
        adata.obs['dca_split'] = 'train'

    adata.obs['dca_split'] = adata.obs['dca_split'].astype('category')
    print('dca: Successfully preprocessed {} genes and {} cells.'.format(adata.n_vars, adata.n_obs))

    return adata

def normalize(adata, filter_min_counts=True, size_factors=True, normalize_input=True, logtrans_input=True):
    if filter_min_counts:
        sc.pp.filter_genes(adata, min_counts=1)
        sc.pp.filter_cells(adata, min_counts=1)

    if size_factors or normalize_input or logtrans_input:
        adata.raw = adata.copy()
    else:
        adata.raw = adata

    if size_factors:
        sc.pp.normalize_per_cell(adata)
        adata.obs['size_factors'] = adata.obs.n_counts / np.median(adata.obs.n_counts)
    else:
        adata.obs['size_factors'] = 1.0

    if logtrans_input:
        sc.pp.log1p(adata)

    if normalize_input:
        sc.pp.scale(adata)

    return adata

  
 
