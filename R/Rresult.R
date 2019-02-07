#' @title Reduced result list
#' @description export reduced isotope labeled result list
#' @param full_Result full result list
#' @import dplyr
#' @return csv file with repeated results being removed
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
#' reduced_result <- Rresult(full_result)

Rresult <- function(full_Result) {
##(1) set global variable
 group_by <- UQ<- filter_at <- vars <- all_vars <- "%>%" <- NULL
 . <- NULL

 ##(2) find base peaks for Labeled molecules
  ## get colnames according to m/z, rt and intensity
  mz_name <- colnames(full_Result)[grepl("mz", colnames(full_Result))]
  rt_name <- colnames(full_Result)[grepl("rt", colnames(full_Result))]
  int_name <- colnames(full_Result)[grepl("intensity", colnames(full_Result))]

  ## search for base peaks
  n_name <- length(mz_name)
  i = 1
  r_result <- vector("list", length = n_name + 1)
  r_result[[1]] <- full_Result
  while (i <= n_name) {
      r_result[[i+1]] <- as.data.frame(r_result[i]) %>%
      group_by(UQ(as.name(mz_name[i])), UQ(as.name(rt_name[i]))) %>%
      filter_at(vars(UQ(int_name[-i])), all_vars(. == max(.)))
      i = i + 1
  }

  ## get reduced peaklist
  reduced_result <- as.data.frame(r_result[[n_name]])
  return(reduced_result)
}

