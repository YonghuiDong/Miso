# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'



# isotope differences, can be extended according to the feeding precusors
H2 = (2.014102 - 1.007825)
C13 = (13.003355 - 12)
N15 = (15.000109 - 14.003074)

# dependency checking, and install required packges if needed
if(!require(S4Vectors)) {
  message("installing the 'S4Vectors' package")
  source("https://bioconductor.org/biocLite.R")
  biocLite("S4Vectors")
}




dual.iso <- function(iso1, n11, n12,  iso2 = 0, n21 = 0, n22 = 0, exp.base,
                     exp.iso, ppm = 100, rt.dif = 0.2) {

  # prepare the labelling pattern according to the feeding precusor. iso2 is 0 by default.
  pattern1 <- seq(n11, n12, by = -1)
  pattern2 <- seq(n21, n22, by = -1)
  iso.pattern1 <- iso1 * pattern1
  iso.pattern2 <- iso2 * pattern2
  pattern <- c(outer(pattern1, pattern2, '+'))
  A <- c(outer(iso.pattern1, iso.pattern2, '+'))
  C <- vector('list', length(A))

  # suppose all the peaks in exp.iso were labeled,
  # then we calculated their correspoinding un-labled peaks
  for (i in 1 : length(A)) {
    C[[i]] <- exp.iso - A[i]
  }

  ## correct RT
  for (i in 1 : length(A)) {
    C[[i]]$RT <- C[[i]]$RT + A[i]
  }

  ## prepare data frame
  C.dframe <- do.call(rbind.data.frame, C)
  z_charge <- rep(pattern, rep(dim(exp.iso)[1],length(A)))
  mz_dif <- rep(A, rep(dim(exp.iso)[1],length(A)))
  C.dframe <- cbind(C.dframe, z_charge, mz_dif)

  # select real isotope labeled peaks according to accepted ppm and RT ranges
  expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL),
                                         list(...))
  CB.expand <- expand.grid.df(exp.base, C.dframe)
  colnames(CB.expand) <-c('B.mz', 'B.rt', 'iso1.mz', 'iso1.rt','charge',
                          'mz_dif')

  ## calculate real ppm
  abs.mz = with(CB.expand, abs(iso1.mz - B.mz) * 10^6 / B.mz)

  ## calculate RT dufference
  abs.rt = with(CB.expand, abs(iso1.rt - B.rt))

  ## select the peaks in our acceprance range
  Result = CB.expand[(abs.mz <= ppm & abs.rt <= rt.dif),]

  ## prepare the result
  Result$iso1.mz = Result$iso1.mz + Result$mz_dif
  Result = Result[, 1:5]
  return(Result)
}



# full result
fResult <- function(iso.C, iso.D) {
  require(S4Vectors)
  p <- as.matrix(findMatches(iso.C$B.mz, iso.D$B.mz))
  fResult <- cbind(Bmz = iso.C$B.mz[p[,1]], Brt = iso.C$B.rt[p[,1]],
                   Cmz = iso.C$iso1.mz[p[,1]], Crt = iso.C$iso1.rt[p[,1]],
                   C_charge = iso.C$charge[p[,1]],
                   Dmz = iso.D$iso1.mz[p[,2]],
                   Drt= iso.D$iso1.rt[p[,2]],
                   D_charge = iso.D$charge[p[,2]])
  return(fResult)
  fResult <- fResult[!duplicated(fResult),]
  write.csv(fResult,'full_result.csv')
}



# simple and dirty result
rResult <- function(full_Result) {
  Result <- as.data.frame(full_Result)
  ## remove B, where both mz and RT are the same
  Result2 <- Result[!duplicated(Result[,c(1,2)]),]
  Result2 <- Result2[order(Result2$Bmz, Result2$Cmz,Result2$Dmz),]
  Result3 <- Result2[!duplicated(Result2[,c(3,4)]),]
  Result4 <- Result3[!duplicated(Result3[,c(6,7)]),]
  write.csv(Result4,'reduced_result.csv')
}


