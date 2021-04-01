# Evaluation procedures for forecasting with spatiotemporal data

This repository contains the research compendium of "Evaluation procedures for forecasting with spatiotemporal data", authored by Mariana Oliveira, Luis Torgo, and Vitor Santos Costa, and published at the Mathematics Journal, on the Special Issue "Spatial Statistics and Its Application".

You are free to use and/or adapt the code we freely provide. However, we do require that if you do that you cite the paper where these results and code were published:

Oliveira M, Torgo L, Santos Costa V. Evaluation Procedures for Forecasting with Spatiotemporal Data. *Mathematics*. 2021; 9(6):691. https://doi.org/10.3390/math9060691

@Article{Oliveira2021, AUTHOR = {Oliveira, Mariana and Torgo, Luís and Santos Costa, Vítor}, TITLE = {Evaluation Procedures for Forecasting with Spatiotemporal Data}, JOURNAL = {Mathematics}, VOLUME = {9}, YEAR = {2021}, NUMBER = {6}, ARTICLE-NUMBER = {691}, URL = {https://www.mdpi.com/2227-7390/9/6/691}, ISSN = {2227-7390}, DOI = {10.3390/math9060691}}


This article is an extended version of a conference paper presented at ECML-PKDD 2018 and published as:

Oliveira M, Torgo L, Santos Costa V. "Evaluation Procedures for Forecasting with Spatio-Temporal Data." In *Proceedings of the European Conference on Machine Learning and Principles and Practice of Knowledge Discovery in Databases, ECML-PKDD* (pp. 703–718). Springer, Cham, 2018. https://doi.org/10.1007/978-3-030-10925-7_43


If you adapt the code to your own needs, you are also required to maintain information on your code concerning the original source of the code (e.g. the URL of this page) and a reference to the original paper.

## Prerequisites

To install this package, run:

```
library(devtools)  # You need to install this package!
install_github("mrfoliveira/STEvaluation-MDPI2021", ref="master")
```

To replicate figures, installing the following package is also necessary:

```
if(!("ggplot2" in installed.packages())) install.packages("ggplot2")
```

## Reproducing experiments

To run experiments, run the following lines from the main directory:

```
library(STEvaluationExt)
PATH <- system.file("inst/", package="STEvaluationExt")

# Artificial experiments
source(paste0(PATH, "/step1_gen_data.R"))
source(paste0(PATH, "/step2_artificial_experiments.R"))

# Real-world experiments
source(paste0(PATH, "/step3_read_data.R"))
source(paste0(PATH, "/step4_get_indicators.R"))
source(paste0(PATH, "/step5_real_experiments.R"))
```