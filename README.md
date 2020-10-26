# famat
Functional analysis of metabolic and transcriptomic data

# Introduction

The aim of `r Rpackage("famat")` is to allow users to determine functional 
links between metabolites and genes. These metabolites and genes lists may 
be related to a specific experiment/study, but `r Rpackage("famat")` only 
needs a gene symbols list and a Kegg Compound ids list. Using these lists, 
`r Rpackage("famat")` performs pathway enrichment analysis, direct interactions
between elements inside pathways extraction, GO terms enrichment analysis, 
calculation of user's elements centrality (number of direct interactions 
between an element and others inside a pathway) and extraction of information
related to user's elements. 

Functions available are: <br><ul>
    <li> path_enrich : pathways enrichment analysis </li>
    <li> interactions : direct interactions and centrality </li>
    <li> compl_data : GO terms enrichment analysis and user's elements data
    extraction </li>
    <li> rshiny : use of previous function's data in a shiny interface </li>
</ul>

# Installation

Run this command line to install `r Rpackage("famat")`.
```{r, eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("famat")
```

Then, load `r Rpackage("famat")` using library.
```{r load_famat, message=FALSE, warning=FALSE}
library(famat)
