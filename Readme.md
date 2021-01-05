# Denoising scRNA-seq Data: A Bayesian Approach

This project explores the use of a deep auto-encoder with a noise model loss function on the task of denoising scRNA-seq data. More information on the method can be found in the [original paper](https://www.nature.com/articles/s41467-018-07931-2). 

# Problem Overview

## A foreword on scRNA-seq Data
Single-cell RNA sequencing (scRNA-seq) is a recent sequencing method capable of getting the counts of expressed genes per cell (as opposed to count per cell sample like [Bulk-seq]() technologies). This new sequencing resolution allows studying biological phenomena at a new resolution. The downside of the method is generating high-dimensional sparse matrices, which are exacerbated by noise-inducing operations like amplification and [add more]. The output can therefore be a large matrix with multiple false zero-counts, known as *dropout*. 

Zero-counts reflect that a gene is repressed for a given cell (i.e., the cell in its current state does not use said gene to function). False zero counts can therefore give the wrong idea about the state of the cell. It would be helpful to identify the false counts and impute them (i.e, extrapolate their values) from neighboring cells and/or genes. This operation is known as *denoising*. Various methods have been proposed for this task, including:
- scImpute:
- 
- And more recently:

This work builds on [DCA]() a method reported to scale linearly with data and improve the results of multiple  scRNA-seq downstream analyses, including cluster separation as shown below.


## Denoising from a Bayesian Perspective
One approach to denoising is to consider the generative model of the counts. It is wildly speculated that read-based methods implicitly assume a [Zero-Inflated Negative Binomial](https://en.wikipedia.org/wiki/Zero-inflated_model)(ZINB) as a count distribution [2]. The generative process for this data looks something like this:

> **input** 𝜋𝑔 ,  𝜇𝑔 ,  𝜃𝑔 ,  𝐺𝑒𝑛𝑒𝑠_{1𝑥𝑔}  // array of genes of size  𝑔   
> **init** 𝑑𝑎𝑡𝑎_{𝑐𝑥𝑔}=0  // matrix of counts of size cells x genes  
> **output** 𝑑𝑎𝑡𝑎_{𝑐𝑥𝑔}  // with updated count values  
> **for** g in  𝐺𝑒𝑛𝑒𝑠_{1𝑥𝑔}:  
> &nbsp;&nbsp;&nbsp;&nbsp; Draw assignments for cells  𝑍𝑐  ~ Bernoulli( 𝜋𝑔 )  
> &nbsp;&nbsp;&nbsp;&nbsp; **for** every cell  𝑐𝑖  in  𝑍𝑐[𝑍𝑐==0] :  
> &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; Set data[ 𝑐𝑖 , g] = 0 (dropout)  
> &nbsp;&nbsp;&nbsp;&nbsp; **for** every cell  𝑐𝑖  in  𝑍𝑐[𝑍𝑐==1] :  
> &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; Set data[ 𝑐𝑖 , g] = mean ( 𝜇𝑔 ) of NB( 𝜇𝑔 ,  𝜃𝑔 )  
> **return 𝑑𝑎𝑡𝑎_{𝑐𝑥𝑔}**  

It can also be visualized through a plate diagram as shown below:

![Plate diagram for DCA](https://github.com/NajwaLaabid/Denoising-sc-RNA-seq/blob/main/img/dca_gm.png "DCA GM")

Using this knowledge, we can fit a [ZINB](https://en.wikipedia.org/wiki/Zero-inflated_model) for every gene in the matrix, and use the mean of the distribution as an imputation value. The fitting (or inference) can be done using a deep auto-encoder model.

## Auto-encoders for inference

# Structure of Repo

```bash
├── data
│   ├── simraw2.pickle
│   ├── simtrue2.pickle
│   ├── simtrue2.pickle
├── img
│   ├── dca_gm.png

```
- report.ipynb: the report of the project. Includes model implementation and experiments.
- original_dca.ipynb: runs the code of the original DCA implementation (available [here](https://github.com/theislab/dca) to be compared with the current (simplified) implementation.
- generate_data.ipynb: notebook to generate simulated data using splatter (available (here)[https://github.com/Oshlack/splatter])
- data/: directory where data files (given data genes.csv and simulated data) are stored. The outputs of the original_dca analysis are also saved there.
- img/: directory holding the images used in the report.
- preprocess.py: helper functions for reading and preprocessing givten data (genes.csv) and performing basic operations like generating plots and normalizing a given anndata object.

# Summary of Results

We tested the model by generating PCA clustering plots on noise-ridden and denoised simulated data. We compared the performance of our model to the original DCA implementation. The results are shown below:


# Reproducibility

We recommend running the notebooks in the following order: generate_data.ipynb > original_dca.ipynb > report.ipynb
The report can also be ran/viewed on its own, since the data generated from the other two notebooks is already available in 'data/'.

# Disclaimer:
This project is based on the Deep Count Auto-encoder (DCA) framework presented in [1]. The code used here is therefore largely adapted from DCA's, with modifications for simplicity and dependencies upgrade (the original DCA code is not compatible with TensorFlow v2 along with other libraries). In particular, the code in generate_data is taken from [here](https://github.com/theislab/dca/blob/master/reproducibility/code/Figure2.ipynb).

# Warnings:
Colab seems to have trouble downloading splatter and its dependencies. We recommend running generate_data.ipynb on a local machine then uploading the generated files to the 'data/' folder in drive if similar issues are faced.

# References

[1] Eraslan, G., Simon, L.M., Mircea, M. et al. Single-cell RNA-seq denoising using a deep count autoencoder. Nat Commun 10, 390 (2019). https://doi.org/10.1038/s41467-018-07931-2

[2] Chen, W., Li, Y., Easton, J. et al. UMI-count modeling and differential expression analysis for single-cell RNA sequencing. Genome Biol 19, 70 (2018). https://doi.org/10.1186/s13059-018-1438-9

