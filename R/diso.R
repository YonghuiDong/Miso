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
#' @param poly polymer of the feeding precursor derived metabolites, e.g. dimer poly = 2, trimer poly =3, default value 1.
#' @return results containing unlabled and their corresponding labled analytes, with RT and labeling information.
#' @export
#' @examples
#' data(lcms)
#' explist <- prefilter(lcms, subgroup = c("B", "C", "D"), unlabel = "B")
#' exp.B <- explist$B
#' exp.C <- explist$C
#' exp.D <- explist$D
#' iso.C <- diso(iso1 = 'H2', n11 = 4, n12 = 2, exp.base = exp.B, exp.iso = exp.C)
#' iso.D <- diso(iso1 = 'C13', n11 = 9, n12 = 6, iso2 = 'N15', n21 = 1, n22 = 0,
#' exp.base = iso.C[,1:3], exp.iso = exp.D)

diso <- function(iso1, n11, n12, iso2 = 'NO', n21 = 0, n22 = 0, exp.base,
                     exp.iso, ppm = 10, rt.dif = 6, poly = 1) {
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
  if (poly <= 0) {stop('poly should be >= 1')}
  if (poly %% 1 != 0) {stop('poly should be integer')}
  cat("done");

  cat("\n(2) Preparing datacube...");
  # prepare the labelling pattern according to the feeding precusor. iso2 is 0 by default.
  pattern1 <- seq(n11, n12, by = -1) * poly
  pattern2 <- seq(n21, n22, by = -1) * poly
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

  ## correct Intensity
  for (i in 1 : length(A)) {
    C[[i]]$intensity <- C[[i]]$intensity + A[i]
  }

  ## prepare data frame
  C.dframe <- do.call(rbind.data.frame, C)
  z_charge <- rep(pattern, rep(dim(exp.iso)[1], length(A)))
  mz_dif <- rep(A, rep(dim(exp.iso)[1],length(A)))
  C.dframe <- cbind(C.dframe, z_charge, mz_dif)
  cat("done");

  cat("\n(3) Performing the 2ed filtering...");
  # select real isotope labeled peaks according to accepted ppm and RT ranges
  expand.grid.df <- function(...) Reduce(function(...) merge(..., by = NULL), list(...))
  CB.expand <- expand.grid.df(exp.base, C.dframe)
  ## name CB.expand according to feeding pattern
  if(iso2 == "NO") {
    iso2 = ""
  } else {
    iso2 = iso2
  }
  colnames(CB.expand) <-c(paste("B.", c("mz", "rt", "intensity"), sep = ""),
                          paste(iso1, iso2, ".", c("mz", "rt", "intensity"), sep = ""),
                          'charge', 'mz_dif')
  ## calculate real ppm
  iso_mz <- paste(iso1, iso2, ".", "mz", sep = "")
  abs.mz = with(CB.expand, abs(CB.expand[iso_mz] - B.mz) * 10^6 / B.mz)

  ## calculate RT dufference
  iso_rt <- paste(iso1, iso2, ".", "rt", sep = "")
  abs.rt = with(CB.expand, abs(CB.expand[iso_rt] - B.rt))

  ## select the peaks in our acceprance range
  Result = CB.expand[(abs.mz <= ppm & abs.rt <= rt.dif),]

  ## prepare the result
  Result[iso_mz] = Result[iso_mz] + Result$mz_dif
  Result$mz_dif = NULL
  cat("done");
  return(Result)
}


