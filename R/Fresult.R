#' @title Full result list
#' @description export full isotope labeled result list
#' @param iso.C result from the first labeled precusor from Exp. C
#' @param iso.D result from the second labled precusor from Exp. D
#' @importFrom S4Vectors findMatches
#' @return file containing the all the possible combined results. Full list, but redundant
#' @export
#' @examples
#' data(lcms)
#' explist <- prefilter(lcms[1: 500, ]) # use a subset of lcms data as example
#' exp.B <- explist$exp.B[, -2]
#' exp.C <- explist$exp.C[, -2]
#' exp.D <- explist$exp.D[, -2]
#' iso.C <- diso(iso1 = 'H2', n11 = 4, n12 = 3, exp.base = exp.B, exp.iso = exp.C)
#' iso.D <- diso(iso1 = 'C13', n11 = 9, n12 = 6, iso2 = 'N15', n21 = 1, n22 = 0,
#' exp.base = iso.C[,1:2], exp.iso = exp.D)
#' full_Result <- Fresult(iso.C, iso.D)
#'
Fresult <- function(iso.C, iso.D) {
  p <- as.matrix(findMatches(iso.C$B.mz, iso.D$B.mz))
  fResult <- cbind(Unlabel_mz = iso.C$B.mz[p[,1]], Unlabel_rt = iso.C$B.rt[p[,1]],
                   Label_I_mz = iso.C$iso1.mz[p[,1]],  Label_I_rt = iso.C$iso1.rt[p[,1]],
                   label_I = iso.C$charge[p[,1]],
                   Label_II_mz = iso.D$iso1.mz[p[,2]],
                   Label_II_rt= iso.D$iso1.rt[p[,2]],
                   Label_II = iso.D$charge[p[,2]])
  fResult <- fResult[!duplicated(fResult),]
  return(fResult)
}
