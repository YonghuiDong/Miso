#' @title Reduced result list
#' @description export reduced isotope labeled result list
#' @param full_Result full result list
#' @param cutint cutoff intensity.Rows that any ion intensities below the cutoff value will removed
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

Rresult <- function(full_Result, cutint = 0) {

  Result <- as.data.frame(full_Result)

  #(1) remove note for "no visible binding for global variable"
  Unlabel_mz <- Unlabel_rt <- Unlabel_int <- Label_I_mz <- Label_I_rt <- Label_I_int <-
    Label_II_mz <- Label_II_rt <- Label_II_int <- NULL

  #(2) find base peaks for Label I and Label II molecules
  Result1 <- Result %>%
    group_by(Unlabel_mz, Unlabel_rt) %>%
    filter(Label_I_int == max(Label_I_int), Label_II_int == max(Label_II_int))

  #(3) find base peaks for unlabled molecules
  ## Filter according to Labeled I
  Result2 <- Result1 %>%
    group_by(Label_I_mz, Label_I_rt) %>%
    filter(Unlabel_int == max(Unlabel_int))
  ## Filter according to Labeled II
  Result3 <- Result2 %>%
    group_by(Label_II_mz, Label_II_rt) %>%
    filter(Unlabel_int == max(Unlabel_int))

  #(4) remove repeated rows that cannot be removed using '!duplicated()' function
  # due to some formatting differences.
  Result4 <- Result3[!duplicated(Result3[,c(1, 2, 3)]),]

  #(5) sort data according to RT
  Result5 <- Result4[order(Result4$Unlabel_rt), ]

  #(6) Filter according to intensity cutoff value. This is used to remove those minor peaks
  # which is difficult to identify
  Result6 <- Result5[Result5$Unlabel_int >= cutint & Result5$Label_I_int >= cutint
                      & Result5$Label_II_int >= cutint, ]
  return(Result6)
}

