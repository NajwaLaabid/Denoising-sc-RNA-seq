Disclaimer:
This project is based on the Deep Count Auto-encoder (DCA) framework presented in [1]. The code used here is therefore largely adapted from DCA's, with modifications for simplicity and dependencies upgrade (the original DCA code is not compatible with TensorFlow v2 along with other libraries). In particular, the code is generate_data is taken from (here)[https://github.com/theislab/dca/blob/master/reproducibility/code/Figure2.ipynb].

Files included:
- report.ipynb: the report of the project. Includes model implementation and experiments.
- original_dca.ipynb: runs the code of the original DCA implementation (available (here)[https://github.com/theislab/dca]) to be compared with the current (simplified) implementation.
- generate_data.ipynb: notebook to generate simulated data using splatter (available (here)[https://github.com/Oshlack/splatter])
- data/: directory where data files (given data genes.csv and simulated data) are stored. The outputs of the original_dca analysis are also saved there.
- img/: directory holding the images used in the report.
- preprocess.py: helper functions for reading and preprocessing givten data (genes.csv) and performing basic operations like generating plots and normalizing a given anndata object.

How to replicate the work:
We recommend running the notebooks in the following order: generate_data.ipynb > original_dca.ipynb > report.ipynb
The report can also be ran/viewed on its own, since the data generated from the other two notebooks is already available in 'data/'.

Warnings:
Colab seems to have trouble downloading splatter and its dependencies. We recommend running generate_data.ipynb on a local machine then uploading the generated files to the 'data/' folder in drive if similar issues are faced.

