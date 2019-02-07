#' @title plot isotopologues
#' @description plot unlabled and labeled Isotopologues from filtering result
#' @param dat isotope filtering result
#' @param rinx row index
#' @import ggplot2
#' @importFrom plotly ggplotly
#' @return interactive plot
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
#' isoplot(full_result, 1)

isoplot <- function(dat, rinx) {
  ##(1) check input
  mz <- int <- NULL
  if(rinx <=0 | rinx > dim(dat)[1])
  {stop(paste('your result index should between 1 and', dim(dat)[1]))}

  ##(2) prepare input
  dat <- as.data.frame(dat)[rinx,]
  rt <- round(dat$B.rt, digits=2)
  mz <- as.numeric(dat[, grepl("mz", colnames(dat))])
  int <- as.numeric(dat[, grepl("intensity", colnames(dat))])
  dat2 <- cbind.data.frame(mz, int)

  ##(3) plot
  p <- ggplot(dat2, aes(x = mz, ymax = int, ymin = 0)) +
    geom_linerange(size = 0.6) +
    xlab("m/z") +
    ylab("Intensity") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(int)*1.1),
                       labels = scales::scientific) +
    theme_bw() +
    ggtitle(paste("Isotopologue at RT = ", rt, "min"))
    ggplotly(p)
}

