#' @title Reduced result list
#' @description export reduced isotope labeled result list
#' @param full_Result full result list
#' @return csv file with repeated results being removed
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
#' reduced_Result <- Rresult(full_Result)

Rresult <- function(full_Result) {
  Result <- as.data.frame(full_Result)
  ## remove B, where both mz and RT are the same
  Result2 <- Result[!duplicated(Result[,c(1,2)]),]
  Result2 <- Result2[order(Result2$Bmz, Result2$Cmz,Result2$Dmz),]
  Result3 <- Result2[!duplicated(Result2[,c(3,4)]),]
  Result4 <- Result3[!duplicated(Result3[,c(6,7)]),]
  return(Result4)
}
