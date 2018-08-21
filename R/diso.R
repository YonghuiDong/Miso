#' @title Isotope filtering
#' @description filtering isotopically labeled analytes according to RT and mass differences
#' @param iso1 the first labled atom in precusor ion.
#' @param n11 the maximum numbers of the first labled atoms expected in the labeled intermediates.
#' @param n12 the minimum numbers of the first labled atoms expected in the labeled intermediates.
#' @param iso2 the second labled atom in the same precusor ion, default value 'NO' (not exist).
#' @param n21 the maximum numbers of the second labled atoms expected in the labeled intermediates, default value 0.
#' @param n22 the minimum numbers of the second labled atoms expected in the labeled intermediates, default value 0.
#' @param exp.base the control group (fed with unlabled precusor).
#' @param exp.iso isotope labeled group.
#' @param ppm m/z tolarance, default value 30.
#' @param rt.dif retention time tolarance, default value 6 seconds.
#' @return results containing unlabled and their corresponding labled analytes, with RT and labeling information.
#' @export
#' @examples
#' data(lcms)
#' explist <- prefilter(lcms[1: 500, ]) # use a subset of lcms data as example
#' exp.B <- explist$exp.B[, -2]
#' exp.C <- explist$exp.C[, -2]
#' exp.D <- explist$exp.D[, -2]
#' iso.C <- diso(iso1 = 'H2', n11 = 4, n12 = 2, exp.base = exp.B, exp.iso = exp.C)
#' iso.D <- diso(iso1 = 'C13', n11 = 9, n12 = 6, iso2 = 'N15', n21 = 1, n22 = 0,
#' exp.base = iso.C[,1:2], exp.iso = exp.D)

diso <- function(iso1, n11, n12, iso2 = 'NO', n21 = 0, n22 = 0, exp.base,
                     exp.iso, ppm = 30, rt.dif = 6) {
  # check the input variables
  cat("\n(1) Checking input parameters...");
  element <- isoDB$element
  if (!is.element(toupper(iso1), element)) {stop('The labeled atom iso1 does not exist.')}
  if (!is.element(toupper(iso2), element)) {stop('The labeled atom iso2 does not exist.')}
  if(!is.numeric(n11)) {stop('n11 is not numeric')}
  if(!is.numeric(n12)) {stop('n12 is not numeric')}
  if(!is.numeric(n21)) {stop('n21 is not numeric')}
  if(!is.numeric(n22)) {stop('n22 is not numeric')}
  if (n11 < n12) {stop('The argument "n11" must be no less than "n12"' )}
  if (n21 < n22) {stop('The argument "n21" must be no less than "n22"' )}
  cat("done");

  cat("\n(2) Preparing datacube...");
  # prepare the labelling pattern according to the feeding precusor. iso2 is 0 by default.
  pattern1 <- seq(n11, n12, by = -1)
  pattern2 <- seq(n21, n22, by = -1)
  mass_dif1 <- isoDB$mass_dif[match(iso1, isoDB$element)]
  mass_dif2 <- isoDB$mass_dif[match(iso2, isoDB$element)]
  iso.pattern1 <- mass_dif1 * pattern1
  iso.pattern2 <- mass_dif2 * pattern2
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
    C[[i]]$rt <- C[[i]]$rt + A[i]
  }

  ## prepare data frame
  C.dframe <- do.call(rbind.data.frame, C)
  z_charge <- rep(pattern, rep(dim(exp.iso)[1],length(A)))
  mz_dif <- rep(A, rep(dim(exp.iso)[1],length(A)))
  C.dframe <- cbind(C.dframe, z_charge, mz_dif)
  cat("done");

  cat("\n(3) Performing the 2ed filtering...");
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
  cat("done");
  return(Result)
}


