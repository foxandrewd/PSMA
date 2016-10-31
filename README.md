# PSMA

R code accompanying the journal article Fox et al. (2016) "Statistical deconvolution enhances interpretability of whole blood DNA methylome data and reveals immune cell type-specific differential methylation in Multiple Sclerosis" 

## Source Files

### psma.R:
This file implements the PSMA statistical deconvolution method for making differential methylation
calls on a cell type specific basis (in Neutrophils, T cells, NK cells, B cells and Monocytes)
from whole-blood-only DNA methylation data

### linear-loocv-R:
This file implements the linear-loocv method for prediction of cellular proportions of
Neutrophils, CD4+ T cells, CD8+ T cells, NK cells, B cells and Monocytes from whole-blood
DNA methylation data

### houseman-code-for-manuscript.R:
This file implements the Houseman deconvolution method for prediction of cellular proportions of
Neutrophils, CD4+ T cells, CD8+ T cells, NK cells, B cells and Monocytes from whole-blood DNA 
methylation data 
