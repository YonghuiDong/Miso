#' @title Prefilter2
#' @description prefiltering isotopically labeled analytes according to the experiment design.
#' @param peak xcms processed dataset.
#' @param nsam number of samples in each experiment, default value = 2. If the sample numbers are different in each
#' gorup, a vector can be used to input the number of samples in each group, i.e. nsam = c(0, 2, 2, 3, 3).
#' @param p p-value threshold, default value = 0.05
#' @param fold fold change threshold, default value = 10
#' @importFrom stats lm
#' @importFrom stats aov
#' @importFrom stats TukeyHSD
#' @return a filtered peaklist
#' @export
#' @examples
#' data(lcms)
#' explist <- prefilter2(lcms[1:100, ])
#'
prefilter2 <- function(peak, nsam = 2, p = 0.05, fold = 10){

  ##(1) check input
  cat("\n(1) Checking input parameters...");
  ## check nsam
  if(!is.numeric(nsam)) {stop("invalid calss of nsam - not numeric")}
  if(length(nsam) != 1 & length(nsam) != 5) {stop("invalid length of nsam, either 1 or 5")}
  if(min(nsam) < 0) {stop("nsam should be no less than 0")}
  if(any(nsam %% 1 != 0)) {stop("nsam should be integer")}
  ## check p value
  if(!is.numeric(p)){stop("invalid calss of p-value threshold: not numeric")}
  if(p > 1 | p < 0){stop("invalid p-value, the value should between 0 and 1")}
  cat("done");

  ##(2) prepare new peaklist
  ## subset peak
  peaklist <- peak[peak$B > 0 | peak$C > 0 | peak$D > 0, ]
  A = peaklist[, c(-1:-12)]
  B = t(A)
  ## add group information
  if(length(nsam) == 1){
    each_exp = rep(nsam,5)
  } else{
    each_exp = nsam
  }
  C = cbind.data.frame(B, group = rep(LETTERS[1:5], each_exp))
  ## only select B, C, D groups
  peaklist_new <- C[(C$group %in% c("B", "C", "D")),]

  ##(3) calculating fold change
  cat("\n(3) Calculating fold change...");
  ret <- fold(peaklist_new[, -dim(peaklist_new)[2]], peaklist_new$group)
  cat("done");

  ##(4) statistical test
  cat("\n(4) Performing statistical test...");
  output <- getp(peaklist_new)
  cat("done");

  ##(5) perform the first filter
  cat("\n(5) Performing the first filtering...");

  ## extract the index in which either p is less than the p-value threshold or fold change
  ## is higher than the fold change threshold
  c_index1 <- which(output$p_CB <= p & output$p_CD <= p)
  c_index2 <- which(ret$F_CB >= fold & ret$F_CD >= fold)
  d_index1 <- which(output$p_BD <= p & output$p_CD <= p)
  d_index2 <- which(ret$F_DB >= fold & ret$F_DC >= fold)
  c_index <- sort(unique(c(c_index1, c_index2)))
  d_index <- sort(unique(c(d_index1, d_index2)))

  peaklistB <- peaklist[peaklist$B > 0, ]
  peaklistC <- peaklist[c_index, ]
  peaklistD <- peaklist[d_index, ]
  exp.B <- cbind.data.frame(mz = peaklistB$mz, intensity = peaklistB[, (13 + each_exp[1])],
                            rt = peaklistB$rt)
  exp.C <- cbind.data.frame(mz = peaklistC$mz, intensity = peaklistC[, (13 + sum(each_exp[1:2]))],
                            rt = peaklistC$rt)
  exp.D <- cbind.data.frame(mz = peaklistD$mz, intensity = peaklistD[, (13 + sum(each_exp[1:3]))],
                            rt = peaklistD$rt)
  stat <- cbind.data.frame(peaklist[, c(1:6)], output, ret)
  exp_list <- list(exp.B = exp.B, exp.C = exp.C, exp.D = exp.D, stat = stat)
  cat("done");
  return(exp_list)
}

