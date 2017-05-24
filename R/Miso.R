#' @title
#' @description
#' @param
#' @importFrom S4Vectors
#' @export


# isotope differences, can be extended according to the feeding precusors
iso_DB <- c(H2 = (2.014102 - 1.007825),
            C13 = (13.003355 - 12),
            N15 = (15.000109 - 14.003074),
            O18 = (17.999160 - 15.994915),
            S34 = (33.967867 - 31.972071),
            0
            )

#' Filtering isotopically labeled analytes according to RT and mass differences
#' @param iso1 the first labled atom in precusor ion
#' @param n11 the maximum numbers of the first labled atoms expected in the labeled intermediates
#' @param n12 the minimum numbers of the first labled atoms expected in the labeled intermediates
#' @param iso2 the second labled atom in the same precusor ion, if exist. it is 0 by default
#' @param n21 the maximum numbers of the second labled atoms expected in the labeled intermediates
#' @param n22 the minimum numbers of the second labled atoms expected in the labeled intermediates
#' @param exp.base the control group (fed with unlabled precusor)
#' @param exp.iso isotope labeled group
#' @param ppm mass tolarance
#' @param rt.dif retention time tolarance
#' @return results containing unlabled and their corresponding labled analytes, with their RT and labeling information

dual.iso <- function(iso1, n11, n12,  iso2 = 0, n21 = 0, n22 = 0, exp.base,
                     exp.iso, ppm = 100, rt.dif = 0.2) {
  # check the input variables
  if (!is.element(iso1, iso_DB)) {
    stop('The labeled atom "iso1" does not exist. The atom letter should be capitalized')
  }

  if (!is.element(iso2, iso_DB)) {
    stop('The labeled atom "iso2" does not exist. The atom letter should be capitalized')
  }

  if (n11 < n12) {
    stop('The argument "n11" must be no less than "n12"' )
  }

  if (n21 < n22) {
    stop('The argument "n21" must be no less than "n22"' )
  }

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



#' Export full result
#' @param iso.C result from the first labeled precusor from Exp. C
#' @param iso.D result from the second labled precusor from Exp. D
#' @return an csv file containing the all the possible combined results. Full list, but redundant
#'
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



#' simple but dirty result
#' @param full_Result the full result
#' @return reduced result
#'
rResult <- function(full_Result) {
  Result <- as.data.frame(full_Result)
  ## remove B, where both mz and RT are the same
  Result2 <- Result[!duplicated(Result[,c(1,2)]),]
  Result2 <- Result2[order(Result2$Bmz, Result2$Cmz,Result2$Dmz),]
  Result3 <- Result2[!duplicated(Result2[,c(3,4)]),]
  Result4 <- Result3[!duplicated(Result3[,c(6,7)]),]
  write.csv(Result4,'reduced_result.csv')
}


