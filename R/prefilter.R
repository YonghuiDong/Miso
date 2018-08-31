#' @title Prefilter
#' @description prefiltering isotopically labeled analytes according to the experiment design.
#' @param peak xcms processed dataset.
#' @param cutint cutoff intensity. Ion intensity below the cutoff value will be considered as noise, default value 0.
#' @param nsam number of samples in each experiment, default value 2. If the sample numbers are different in each
#' gorup, a vector can be used to input the number of samples in each group, i.e. nsam = c(0, 2, 2, 3, 3).
#' @param minsam minimum number of samples.The peak is considered valid only when it is at least detected in
#' the minimum number of samples in a group, default value 1.  A vector can be used to To control the minsam in each
#' group, i.e. nsam = c(1, 1, 2, 2, 1).
#' @return a filtered peaklist.
#' @export
#' @examples
#' data(lcms)
#' explist <- prefilter(lcms[1:100, ])

prefilter <- function(peak, cutint = 0, nsam = 2, minsam = 1){
  cat("\n(1) Checking input parameters...");
  # check cutint
  if(!is.numeric(cutint)) {stop("invalid class of cutint - not numeric")}
  if(cutint < 0) {stop("cut off intensity value should be no less than 0")}
  # check nsam
  if(!is.numeric(nsam)) {stop("invalid calss of nsam - not numeric")}
  if(length(nsam) != 1 & length(nsam) != 5) {stop("invalid length of nsam, either 1 or 5")}
  if(min(nsam) < 0) {stop("nsam should be no less than 0")}
  if(any(nsam %% 1 != 0)) {stop("nsam should be integer")}
  # check minsam
  if(!is.numeric(minsam)){stop("invalid calss of minsam - not numeric")}
  if(length(minsam) != 1 & length(minsam) != 5){stop("invalid length of minsam, either 1 or 5")}
  if(min(minsam) < 0){stop("nsam should be no less than 0")}
  if(any(minsam %% 1 != 0)){stop("minsam should be integer")}
  # make sure minsam is no over than nsam
  # experiment number in each experiment group
  if(length(nsam) == 1){
    each_exp = rep(nsam,5)
  } else{
    each_exp = nsam
  }
  if(any(nsam >= minsam) == F){stop("minsam should be no over than nsam")}
  cat("done");

  cat("\n(2) Prepareing datacube...");
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
  cat("done");

  ## calculating fold change
  cat("\n(3) Calculating fold change...");
  ret <- fold(peaklist_new[, -dim(peaklist_new)[2]], peaklist_new$group)
  cat("done");

  cat("\n(4) Performing the first filtering...");
  # first filtering & prepare datacube. Regarding intensity, only the intensity of the first sample in each Exp is used
  if (cutint == 0){
    # calculate the fitted index from each group
    b_index <- which(peaklist$B >= minsam)
    c_index <- which(peaklist$C >= minsam & peaklist$D == 0 & peaklist$B == 0)
    d_index <- which(peaklist$D >= minsam & peaklist$C == 0 & peaklist$B == 0)
    # get fitted peaklist
    peaklistB <- peaklist[b_index, ]
    peaklistC <- peaklist[c_index, ]
    peaklistD <- peaklist[d_index, ]
    # add mz, int and rt
    exp.B <- cbind.data.frame(mz = peaklistB$mz, intensity = ret[b_index, ]$B, rt = peaklistB$rt)
    exp.C <- cbind.data.frame(mz = peaklistC$mz, intensity = ret[c_index, ]$C, rt = peaklistC$rt)
    exp.D <- cbind.data.frame(mz = peaklistD$mz, intensity = ret[d_index, ]$D, rt = peaklistD$rt)
    stat <- cbind.data.frame(peaklist[, c(1:6)], ret)
    exp_list <- list(exp.B = exp.B, exp.C = exp.C, exp.D = exp.D, stat = stat)
  } else{
    m = sum(each_exp)
    n = dim(peaklist)[2]
    # thoese intensity values less than cutint will be set to 0
    peaklist[, (n-m+1):n][peaklist[, (n-m+1):n] < cutint] <- 0
    # need to recalculate B, C, D, E
    # creat empty vectors
    B <- vector('integer', length = dim(peaklist)[1])
    C <- vector('integer', length = dim(peaklist)[1])
    D <- vector('integer', length = dim(peaklist)[1])

    # calculate the beginning and end intensity index for each Exp.
    # here considers the sample sizes among each Exp can be different
    b1 <- 13 + each_exp[1]
    b2 <- 12 + sum(each_exp[1:2])
    c1 <- 13 + sum(each_exp[1:2])
    c2 <- 12 + sum(each_exp[1:3])
    d1 <- 13 + sum(each_exp[1:3])
    d2 <- 12 + sum(each_exp[1:4])
    # The intensity over than 0 is transformed into 1
    peaklist2 <-peaklist
    peaklist2[, (n-m+1):n][peaklist2[, (n-m+1):n] > 0] <- 1
    for (i in 1: dim(peaklist2)[1]){
      # calculate new B, C, D, E
      B[i] <- sum(peaklist2[, b1][i] : peaklist2[, b2][i])
      C[i] <- sum(peaklist2[, c1][i] : peaklist2[, c2][i])
      D[i] <- sum(peaklist2[, d1][i] : peaklist2[, d2][i])
    }
    # replace old peaklist B, C, D, E
    peaklist$B <- B
    peaklist$C <- C
    peaklist$D <- D
    # calculate the fitted index from each group
    b_index <- which(peaklist$B >= minsam)
    c_index <- which(peaklist$C >= minsam & peaklist$D == 0 & peaklist$B == 0)
    d_index <- which(peaklist$D >= minsam & peaklist$C == 0 & peaklist$B == 0)
    # get fitted peaklist
    peaklistB <- peaklist[b_index, ]
    peaklistC <- peaklist[c_index, ]
    peaklistD <- peaklist[d_index, ]
    # add mz, int and rt
    exp.B <- cbind.data.frame(mz = peaklistB$mz, intensity = ret[b_index, ]$B, rt = peaklistB$rt)
    exp.C <- cbind.data.frame(mz = peaklistC$mz, intensity = ret[c_index, ]$C, rt = peaklistC$rt)
    exp.D <- cbind.data.frame(mz = peaklistD$mz, intensity = ret[d_index, ]$D, rt = peaklistD$rt)
    stat <- cbind.data.frame(peaklist[, c(1:6)], ret)
    exp_list <- list(exp.B = exp.B, exp.C = exp.C, exp.D = exp.D, stat = stat)
  }
  cat("done");
  return(exp_list)
}
