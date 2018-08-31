#' @title Reduced result list
#' @description export reduced isotope labeled result list
#' @param full_Result full result list
#' @return csv file with repeated results being removed
#' @export
#' @examples
#' data(lcms)
#' explist <- prefilter(lcms[1: 100, ]) # use a subset of lcms data as example
#' exp.B <- explist$exp.B
#' exp.C <- explist$exp.C
#' exp.D <- explist$exp.D
#' iso.C <- diso(iso1 = 'H2', n11 = 4, n12 = 3, exp.base = exp.B, exp.iso = exp.C)
#' iso.D <- diso(iso1 = 'C13', n11 = 9, n12 = 6, iso2 = 'N15', n21 = 1, n22 = 0,
#' exp.base = iso.C[,1:3], exp.iso = exp.D)
#' full_Result <- Fresult(iso.C, iso.D)
#' reduced_Result <- Rresult(full_Result)

Rresult <- function(full_Result) {
  Result <- as.data.frame(full_Result)
  ## remove B, where mz, Intensity, and RT are the same
  Result2 <- Result[!duplicated(Result[,c(1, 2, 3)]),]
  Result2 <- Result2[order(Result2$Unlabel_mz, Result2$Label_I_mz,Result2$Label_II_mz),]
  Result3 <- Result2[!duplicated(Result2[,c(4, 5, 6)]),]
  Result4 <- Result3[!duplicated(Result3[,c(8, 9, 10)]),]
  return(Result4)
}
