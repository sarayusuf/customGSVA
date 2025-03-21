# customGSVA

An R package extending GSVA with AUCell integration for single-sample gene set enrichment analysis.

## Installation
```R
devtools::install_github("sarayusuf/customGSVA")

library(customGSVA)
data(gset.idx.list)
expr <- matrix(rnorm(1000), nrow = 100, ncol = 10)

# Run GSVA
es_gsva <- gsva(expr, gset.idx.list, method = "gsva")

# Run AUCell
es_aucell <- gsva(expr, gset.idx.list, method = "aucell")

