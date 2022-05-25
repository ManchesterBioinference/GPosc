# GPosc

This is the GitHub page for the article:

Identifying stochastic oscillations in single-cell live imaging time series using Gaussian processes 

Citation: Phillips NE, Manning C, Papalopulu N, Rattray M (2017) Identifying stochastic oscillations in single-cell live imaging time series using Gaussian processes. PLoS Comput Biol 13(5): e1005479. https://doi.org/10.1371/journal.pcbi.1005479

The original article used MATLAB, but there is a reimplementation using Python and GPflow.

## Python

Inside /GPosc/Python there is a Jupyter Notebook tutorial (tutorial.ipynb) that describes the main steps of the pipeline.

There is also a Python script (PaperResultsFig8.py) that generates the main results of the paper (the number of Hes1 or control promoter cells that are classified as oscillating) and Figure 8.

The Conda environment required to run these two scripts is saved as /GPosc/Python/GPosc.yml

The data required to run these files is found in the /GPosc/Python/Data folder.

## Matlab

These are the original MATLAB files used to create the figures in the paper.

The package is based on GPML, and the README for the GPML package is included as READMEGPML. The GPML files are zipped into the file "GPML files.zip". These need to be placed into your working MATLAB folder.

Note that it is important to run the file "startup" before using GPML.

The reference for the GPML toolbox is "C Rasmussen and N Hannes. Gaussian processes for machine learning (GPML) toolbox. Journal of Machine Learning Research, 11:3011 3015, 2010."
