#' @title Prefilter
#' @description prefiltering isotopically labeled analytes according to the experiment design.
#' @param peaklist xcms processed dataset.
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
#' explist <- prefilter(lcms)

prefilter <- function(peaklist, cutint = 0, nsam = 2, minsam = 1){
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
  if(min(minsam) <= 0){stop("nsam should be over than 0")}
  if(any(minsam %% 1 != 0)){stop("minsam should be integer")}
  # make sure minsam is no over than nsam
  # experiment number in each experiment group
  if(length(nsam) == 1){
    each_exp = rep(nsam,5)
  } else{
    each_exp = nsam
  }
  # min number in each experiment group
  if(length(minsam) == 1){
    min_sam = rep(minsam,5)
  } else{
    min_sam = minsam
  }
  if(any(nsam >= minsam) == F){stop("minsam should be no over than nsam")}
  cat("done.");
  cat("\n(2) Performing the first filtering...");
  # experiment number in each experiment group
  if(length(nsam) == 1){
    each_exp = rep(nsam,5)
  } else{
    each_exp = nsam
  }
  # min number in each experiment group
  if(length(minsam) == 1){
    min_sam = rep(minsam,5)
  } else{
    min_sam = minsam
  }
  # replce NA intensity values with 0 from xcms peaklist
  m = sum(each_exp)
  n = dim(peaklist)[2]
  peaklist[, (n-m+1):n][is.na(peaklist[, (n-m+1):n])] <- 0
  # first filtering & prepare datacube. Regarding intensity, only the intensity of the first sample in each Exp is used
  if (cutint == 0){
    peaklistB <-peaklist[peaklist$B >= min_sam[2],]
    peaklistC <- peaklist[peaklist$E >= min_sam[5] & peaklist$C >= min_sam[3] & peaklist$D == 0 & peaklist$B == 0 ,]
    peaklistD <- peaklist[peaklist$E >= min_sam[5] & peaklist$D >= min_sam[4] & peaklist$C == 0 & peaklist$B == 0 ,]
    exp.B <- cbind.data.frame(mz = peaklistB$mz, intensity = peaklistB[, (13 + each_exp[1])],
                              rt = peaklistB$rt)
    exp.C <- cbind.data.frame(mz = peaklistC$mz, intensity = peaklistC[, (13 + sum(each_exp[1:2]))],
                              rt = peaklistC$rt)
    exp.D <- cbind.data.frame(mz = peaklistD$mz, intensity = peaklistD[, (13 + sum(each_exp[1:3]))],
                              rt = peaklistD$rt)

    exp_list <- list(exp.B = exp.B, exp.C = exp.C, exp.D = exp.D)
  } else{
    # thoese intensity values less than cutint will be set to 0
    peaklist[, (n-m+1):n][peaklist[, (n-m+1):n] < cutint] <- 0
    # need to recalculate B, C, D, E
    # creat empty vectors
    B <- vector('integer', length = dim(peaklist)[1])
    C <- vector('integer', length = dim(peaklist)[1])
    D <- vector('integer', length = dim(peaklist)[1])
    E <- vector('integer', length = dim(peaklist)[1])
    # calculate the beginning and end intensity index for each Exp.
    # here considers the sample sizes among each Exp can be different
    b1 <- 13 + each_exp[1]
    b2 <- 12 + sum(each_exp[1:2])
    c1 <- 13 + sum(each_exp[1:2])
    c2 <- 12 + sum(each_exp[1:3])
    d1 <- 13 + sum(each_exp[1:3])
    d2 <- 12 + sum(each_exp[1:4])
    e1 <- 13 + sum(each_exp[1:4])
    e2 <- 12 + sum(each_exp[1:5])
    # The intensity over than 0 is transformed into 1
    peaklist2 <-peaklist
    peaklist2[, (n-m+1):n][peaklist2[, (n-m+1):n] > 0] <- 1
    for (i in 1: dim(peaklist2)[1]){
      # calculate new B, C, D, E
      B[i] <- sum(peaklist2[, b1][i] : peaklist2[, b2][i])
      C[i] <- sum(peaklist2[, c1][i] : peaklist2[, c2][i])
      D[i] <- sum(peaklist2[, d1][i] : peaklist2[, d2][i])
      E[i] <- sum(peaklist2[, e1][i] : peaklist2[, e2][i])
    }

    # replace old peaklist B, C, D, E
    peaklist$B <- B
    peaklist$C <- C
    peaklist$D <- D
    peaklist$E <- E
    # perform the prefilter with cutint is not equal to 0
    peaklistB <-peaklist[peaklist$B >= min_sam[2],]
    peaklistC <- peaklist[peaklist$E >= min_sam[5] & peaklist$C >= min_sam[3] & peaklist$D == 0 & peaklist$B == 0 ,]
    peaklistD <- peaklist[peaklist$E >= min_sam[5] & peaklist$D >= min_sam[4] & peaklist$C == 0 & peaklist$B == 0 ,]
    exp.B <- cbind.data.frame(mz = peaklistB$mz, intensity = peaklistB[, (13 + each_exp[1])],
                              rt = peaklistB$rt)
    exp.C <- cbind.data.frame(mz = peaklistC$mz, intensity = peaklistC[, (13 + sum(each_exp[1:2]))],
                              rt = peaklistC$rt)
    exp.D <- cbind.data.frame(mz = peaklistD$mz, intensity = peaklistD[, (13 + sum(each_exp[1:3]))],
                              rt = peaklistD$rt)

    exp_list <- list(exp.B = exp.B, exp.C = exp.C, exp.D = exp.D)
  }
  cat("done.");
  return(exp_list)
}
