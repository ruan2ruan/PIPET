# PIPET

Phenotypic information based on bulk data predicts relevant subpopulations in single cell data

## Introduction

`PIPET` can be used to predict relevant subpopulations in single-cell data from phenotypic information in bulk data. You can use known feature vectors of phenotypic information to preform PIPET() or create feature vectors via the function in the `PIPET`. This package also provides commonly used downstream analysis for creating customizable visualization results. The workflow of `PIPET` is shown in the following Figure:

<p align="center">
<img src=Figure_PIPET.jpg height="900" width="640">
</p>

## Installation

To install this package, start R (*version >= 4.3.1*) and enter:

``` {r}
if (!require("devtools")) 
    install.packages("devtools")
devtools::install_github("ruan2ruan/PIPET")
```

## Tutorial

In the R terminal, please use the command `?PIPET` to access the help documents.

To view the detailed guide of how to use PIPET, please find package guidance in [vignette](https://ruan2ruan.github.io/PIPET.html), or run the following lines in R:

```{r}
library(PIPET)
browseVignettes("PIPET")
```

## Issues

Any questions about PIPET could be posted to the Issue section of GitHub homepage at https://github.com/ruan2ruan/PIPET/issues.

