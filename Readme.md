# Denoising scRNA-seq Data: A Bayesian Approach

This project explores the use of a deep auto-encoder with a noise model loss function on the task of denoising scRNA-seq data. More information on the method can be found in the [original paper](https://www.nature.com/articles/s41467-018-07931-2). 

# Problem Overview

## A foreword on scRNA-seq Data
Single-cell RNA sequencing (scRNA-seq) is a recent sequencing method capable of getting the counts of expressed genes per cell (as opposed to count per cell sample like [Bulk-seq]() technologies). This new sequencing resolution allows studying biological phenomena at a new resolution. The downside of the method is generating high-dimensional sparse matrices, which are exacerbated by noise-inducing operations like amplification and [add more]. The output can therefore be a large matrix with multiple false zero-counts, known as *dropout*. 

Zero-counts reflect that a gene is repressed for a given cell (i.e., the cell in its current state does not use said gene to function). False zero counts can therefore give the wrong idea about the state of the cell. It would be helpful to identify the false counts and impute them (i.e, extrapolate their values) from neighboring cells and/or genes. This operation is known as *denoising*. 

This work builds on [DCA](https://github.com/theislab/dca/) a method reported to scale linearly with data and improve the results of multiple  scRNA-seq downstream analyses, including cluster separation as shown below.


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

A deep auto-encoder model is proposed as a mechanism for infering the parameters of the noise model. The model is shown to scale linearly with the size of the data [1], which is an important advantage especially when dealing with shallow scRNA-seq data. The model used in this study has a fixed architecture. More extensive hyper parameter/architecture tuning (and flexibility in setting both) as proposed in the [original work](https://github.com/theislab/dca/).

# Structure of Repo

* `data/`: holds the simulated data for experiments. 
* `img/`: directory holding the images used in the report.
* `model/`: code for building the AE model. The actual model is defined in `zae.py`.
* `utils/`: helper code for reading the data and generating plots.
	* `data_utils.py`: helper code for handling the data.
	* `plot_utils.py`: helper functions for plots.
	* `simulate.r`: an *R* script to simulate scRNA-seq data using [splatter](https://github.com/Oshlack/splatter) library. Code is taken from an [example](https://github.com/theislab/dca/blob/master/reproducibility/code/Figure2.ipynb) notebook in the dca repo.
* `experiments.ipynb`: shows how denoising helps recover the original clusters of a noise-ridden simulated data.
* `original_dca.ipynb`: runs the code of the original dca library.

# Summary of Results

We tested the model by generating PCA clustering plots on noise-ridden and denoised simulated data. We compared the performance of our model to the original DCA implementation. The results are shown below:

![Denoising results](https://github.com/NajwaLaabid/Denoising-sc-RNA-seq/blob/main/img/compare_dca.png "Results")

# Warnings & Future Work
* A GPU is recommended for training the model.
* The [DCA](https://github.com/theislab/dca/) library has old dependencies. The notebook executing DCA code installs the specific versions required of the libraries in the first line.
* TODO: figure out why Google Colab seems to have a bit of trouble processing the *R* simulation script. It is recommended to run generate data on your local machine for now.
* TODO: rewrite the loss function to be compatible with TF 2.0. Currently disabling eager execution provides a quick fix to the problem (as shown in `experiments.ipynb`).

# References

[1] Eraslan, G., Simon, L.M., Mircea, M. et al. Single-cell RNA-seq denoising using a deep count autoencoder. Nat Commun 10, 390 (2019). https://doi.org/10.1038/s41467-018-07931-2

[2] Chen, W., Li, Y., Easton, J. et al. UMI-count modeling and differential expression analysis for single-cell RNA sequencing. Genome Biol 19, 70 (2018). https://doi.org/10.1186/s13059-018-1438-9

