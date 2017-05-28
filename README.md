Miso

## Description

- A fast and efficient approach for fishing out the dual isotope labeled analytes.
- Can be easily extended to multiple isotope labeling experiments

## Attention

This packgage is xcms based. xcms processing is crucial for the subsequent processing, either false positive or false negative results from xcms can strongly affect the number of selected labeled analytes. It is therefore strongly recommended to optimize the xcms parameters.

## Usage

> devtools::install_github("YonghuiDong/Miso")
rm(list = ls())
Library(Miso)
mydata <- read.csv('neg_dong_parameter.csv', header = T)
mydata$rt <- mydata$rt/60 # convert second to min

# 1st filtering, according to experiment design

mydataB <-mydata[mydata$B != 0,]

mydataC <- mydata[mydata$E !=0 & mydata$C != 0 & mydata$D == 0 
                  & mydata$B == 0 ,]
mydataD <- mydata[mydata$E !=0 & mydata$D != 0 & mydata$C == 0 
                  & mydata$B == 0 ,]

exp.B <- cbind.data.frame(mz = mydataB$mz, RT = mydataB$rt)
exp.C <- cbind.data.frame(mz = mydataC$mz, RT = mydataC$rt)
exp.D <- cbind.data.frame(mz = mydataD$mz, RT = mydataD$rt)

# group C was fed with H
iso.C <- dual.iso(iso1 = H2, n11 = 5, n12 = 3, exp.base = exp.B, 
                  exp.iso = exp.C)

# Group D was fed with C13, and N15, here we are in
iso.D <- dual.iso(iso1 = C13, n11 = 9, n12 = 6, iso2 = N15, n21 = 1, n22 = 0,
                  exp.base = iso.C[,1:2], exp.iso = exp.D)

full_Result <- fResult(iso.C, iso.D)
reduced_Result <- rResult(full_Result)

