# SiNMFiD

SiNMFiD (Supervised iNMF-informed Deconvolution) is a package implementing the algorithm of the same name (described in [Mesoscale Properties of Molecular Cell Types in the Mouse Brain" by Kriebel et al. (2024)](https://docs.google.com/document/d/166X4o_6HegeS0uUHRa4ahKedkaVEMeiRagCm7T9nOVU/edit?usp=sharing)) for cell-type deconvolution of 3D spatial transcriptomic data with multiple single cell sequencing references. The package also provides utilities for managing multiple parallel analyses for tissue atlas generation and for statistics and plotting.

## Installation

The package is developed and tested under R>=4.3.0. Users can install R following the [instruction provided on CRAN](https://cran.r-project.org/). [RStudio](https://posit.co/downloads/) is a recommended IDE for working with R projects. 

To install SiNMFiD in R, run the following command in an R console:

```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("jsodicoff/SiNMFiD")
```

## Vignette

We have provided a vignette describing how to reproduce and iterate on the results described in "Mesoscale Properties of Molecular Cell Types in the Mouse Brain" by Kriebel et al. (2024) [here](https://github.com/jsodicoff/SiNMFiD/blob/main/vignettes/vignette.html).  
