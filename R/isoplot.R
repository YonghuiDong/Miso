#' @title plot isotopologues
#' @description plot unlabled and labeled Isotopologues from filtering result
#' @param dat isoto filtering result
#' @param rinx row index
#' @import plotly
#' @return interactive plot
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
#' isoplot(full_Result, 1)

isoplot <- function(dat, rinx) {
  ##(1) check input
  mz <- int <- NULL
  if(rinx <=0 | rinx > dim(dat)[1])
  {stop(paste('your result index should between 1 and', dim(dat)[1]))}

  ##(2) prepare input
  dat <- as.data.frame(dat)[rinx,]
  dat2 <- cbind.data.frame(mz = c(dat$Unlabel_mz, dat$Label_I_mz, dat$Label_II_mz),
                           Int = c(dat$Unlabel_int, dat$Label_I_int, dat$Label_II_int))
  dat2[is.na(dat2)] = 0

  ##(3) plot
  plot_ly() %>%
    add_lines(x = rep(dat$Unlabel_mz, 2), y = c(0, dat$Unlabel_int), name = "unlabeled",
              line = list(color = "black", width = 2)) %>%
    add_lines(x = rep(dat$Label_I_mz, 2), y = c(0, dat$Label_I_int), name = "labeled I",
              line = list(color = "red", width = 2)) %>%
    add_lines(x = rep(dat$Label_II_mz, 2), y = c(0, dat$Label_II_int), name = "labeled II",
              line = list(color = "blue", width = 2)) %>%
    layout(
      title = paste("Isotopologues @ RT=", dat$Unlabel_rt, "min"),
      xaxis = list(
        title = 'm/z',
        range = c(min(dat2$mz) - 3, max(dat2$mz) + 3)
      ),
      yaxis = list(
        title = 'Mean Intensity',
        range = c(0, max(dat2$Int) * 1.1)
      )
    )
}
