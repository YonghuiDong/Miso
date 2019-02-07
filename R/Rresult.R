#' @title Reduced result list
#' @description export reduced isotope labeled result list
#' @param full_Result full result list
#' @return csv file with repeated results being removed
#' @export
#' @examples
#' a = 1+1

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

