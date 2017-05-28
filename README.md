Miso: Multi-isotope Labeling for Metabolomics Analysis

## Description

- A fast and efficient approach for fishing out the dual isotope labeled analytes
- Can be easily extended to multiple isotope labeling experiments


## Usage

Example 

```r 
devtools::install_github("YonghuiDong/Miso") 
rm(list = ls())
Library(Miso)
mydata <- read.csv('neg_dong_parameter.csv', header = T)
mydata$rt <- mydata$rt/60 ## convert second to min

## 1st filtering, according to the experiment design
mydataB <-mydata[mydata$B != 0,]
mydataC <- mydata[mydata$E !=0 & mydata$C != 0 & mydata$D == 0 
                  & mydata$B == 0 ,]
mydataD <- mydata[mydata$E !=0 & mydata$D != 0 & mydata$C == 0 
                  & mydata$B == 0 ,]

## Prepare the Datacube
exp.B <- cbind.data.frame(mz = mydataB$mz, RT = mydataB$rt)
exp.C <- cbind.data.frame(mz = mydataC$mz, RT = mydataC$rt)
exp.D <- cbind.data.frame(mz = mydataD$mz, RT = mydataD$rt)

## 2ed filtering, according to the labeling patterns of interest

## group C was fed with H
iso.C <- dual.iso(iso1 = H2, n11 = 5, n12 = 3, exp.base = exp.B, 
                  exp.iso = exp.C)

## Group D was fed with C13, and N15
iso.D <- dual.iso(iso1 = C13, n11 = 9, n12 = 6, iso2 = N15, n21 = 1, n22 = 0,
                  exp.base = iso.C[,1:2], exp.iso = exp.D)

## Generate results
full_Result <- fResult(iso.C, iso.D)
reduced_Result <- rResult(full_Result)
```
## Attention

1. This packgage is `xcms` based. Either false positive or false negative results from xcms can bias the number of selected labeled analytes. It is strongly recommended to optimize the xcms parameters.     

2. R memory limit error might appear during data processing especially for high resolution dataset:   

`Error: memory exhausted (limit reached?), Error during wrapup: memory exhausted (limit reached?)` 

This error is due to the following script:

```r
## In group C, we are looking for analytes labeled with 5, 4, or 3 deuterium (H2).
iso.C <- dual.iso(iso1 = H2, n11 = 5, n12 = 3, exp.base = exp.B, 
                  exp.iso = exp.C)
```

To solve this memory limit problem, this script can be decomposed into 3 steps, which respectively search for analytes labled with 5, 4, or 3 deuterium (H2).

```r
iso.C5 <- dual.iso(iso1 = H2, n11 = 5, n12 = 5, exp.base = exp.B, 
                  exp.iso = exp.C)
iso.C4 <- dual.iso(iso1 = H2, n11 = 4, n12 = 4, exp.base = exp.B, 
                  exp.iso = exp.C)
iso.C3 <- dual.iso(iso1 = H2, n11 = 3, n12 = 3, exp.base = exp.B, 
                  exp.iso = exp.C)

## The results are then combined as iso.C:
iso.C <- rbind(iso.C5, iso.C4, iso.C3)
```

The decomposed step is usually only necessasy for iso.C, as the result list has been significantly reduced. we do not have to do it again for iso.D.
