import random
import anndata
import numpy as np
import tensorflow as tf
import os
import scanpy as sc

from utils import data_utils
from model import zae

def denoise(adata):
    assert isinstance(adata, anndata.AnnData), 'adata must be an AnnData instance'

    # set seed for reproducibility
    random_state=0
    random.seed(random_state)
    np.random.seed(random_state)
    tf.random.set_seed(random_state)
    os.environ['PYTHONHASHSEED'] = '0'

    # data manipulation
    # this creates adata.raw with raw counts and copies adata if copy==True
    adata = data_utils.read_dataset(adata, test_split=True, copy=True)

    # check for zero genes => why not handle in normalize
    nonzero_genes, _ = sc.pp.filter_genes(adata.X, min_counts=1)
    assert nonzero_genes.all(), 'Please remove all-zero genes before using DCA.'

    adata = data_utils.normalize(adata,
                      filter_min_counts=False, # no filtering, keep cell and gene idxs same
                      size_factors=True,
                      normalize_input=True,
                      logtrans_input=False)
    
    network_kwds = {
        'input_size': adata.n_vars,
        'output_size': adata.n_vars,
        'hidden_size': (64, 32, 64),
        'hidden_dropout': 0.,
        'batchnorm': True,
        'activation': 'relu',
        'init': 'glorot_uniform'
    }

    net = zae.ZINBAutoencoder(**network_kwds)
    net.save()
    net.build()

    training_kwds = {
        'epochs': 300,
        'reduce_lr': 10,
        'early_stop': 15,
        'batch_size': 32,
        'optimizer': 'rmsprop',
        'verbose': False,
        'threads': 1,
    }

    hist = net.train(adata[adata.obs.dca_split == 'train'], **training_kwds)

    res = net.predict(adata, return_info=True, copy=True)

    copy=True
    adata = res if copy else adata

    return_info=True
    if return_info:
        adata.uns['dca_loss_history'] = hist.history

    return_model=False
    if return_model:
        return (adata, net) if copy else net
    else:
        return adata if copy else None