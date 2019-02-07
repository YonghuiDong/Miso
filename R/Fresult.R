#' @title Full result list
#' @description export full isotope labeled result list
#' @param x results generated from diso() function. The input can be multiple depending on the
#' experiment design: i.e. single labeling, dual labeling and multiple lableing.
#' @param ... results generated from diso() function
#' @return file containing the all the possible combined results.
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
#' full_result <- Fresult(iso.C, iso.D)

Fresult <- function(x, ...) {
  ## merge result
  p <- Reduce(function(x, y) merge(x, y, by = c("B.mz", "B.rt", "B.intensity")), list(x, ...))
  p <- p[!duplicated(p), ]
  rt_inx <- grepl("rt", colnames(p))
  p[rt_inx] <- round(p[rt_inx]/60, digits=2)
  return(p)
}
