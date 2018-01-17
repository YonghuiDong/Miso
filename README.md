Miso: Multi-Isotope Labeling for Metabolomics Analysis

Version: 0.1.1

## Description

- An efficient approach for fishing out the dual or multiple isotope labeling assisted metabolomics data analytes


## Usage

Example 

```r 
install.packages("Miso")
Library(Miso)
data(lcms)
## First filtering, according to the experiment design
explist <- prefilter(lcms)
## Second filtering, according to the labeling patterns
## group C was fed with H2
iso.C <- diso(iso1 = 'H2', n11 = 5, n12 = 3, exp.base = exp.B, 
                  exp.iso = exp.C)
## Group D was fed with C13, and N15
iso.D <- diso(iso1 = 'C13', n11 = 9, n12 = 6, iso2 = 'N15', n21 = 1, n22 = 0,
                  exp.base = iso.C[,1:2], exp.iso = exp.D)
## Generate results
full_Result <- Fresult(iso.C, iso.D)
reduced_Result <- Rresult(full_Result)
```
## Attention    

1. R memory limit error may appear during data processing especially for high resolution dataset:   

`Error: memory exhausted (limit reached?), Error during wrapup: memory exhausted (limit reached?)` 

This error is due to the following script:

```r
## In group C, we are looking for analytes labeled with 5, 4, or 3 deuteriums (H2).
iso.C <- diso(iso1 = H2, n11 = 4, n12 = 2, exp.base = exp.B, 
                  exp.iso = exp.C)
```

To solve this memory limit problem, the above script can be decomposed into 3 sub-scripts, which respectively search for analytes labled with 4, 3, and 2 deuteriums (H2).

```r
iso.C5 <- diso(iso1 = 'H2', n11 = 4, n12 = 4, exp.base = exp.B, 
                  exp.iso = exp.C)
iso.C4 <- diso(iso1 = 'H2', n11 = 3, n12 = 3, exp.base = exp.B, 
                  exp.iso = exp.C)
iso.C3 <- diso(iso1 = 'H2', n11 = 2, n12 = 2, exp.base = exp.B, 
                  exp.iso = exp.C)

## The results are then combined as iso.C:
iso.C <- rbind(iso.C5, iso.C4, iso.C3)
```

The decomposition step is only usually necessasy for iso.C, as the result list has been significantly reduced. we do not have to do it again for iso.D.
