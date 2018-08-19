#' @title Prefilter2
#' @description prefiltering isotopically labeled analytes according to the experiment design.
#' @param peaklist xcms processed dataset.
#' @param nsam number of samples in each experiment, default value = 2. If the sample numbers are different in each
#' gorup, a vector can be used to input the number of samples in each group, i.e. nsam = c(0, 2, 2, 3, 3).
#' @param p p-value threshold, default value = 0.05
#' @return a filtered peaklist.
#' @export
#' @examples
#' data(lcms)[1:100, ]
#' explist <- prefilter2(lcms)
#' 
prefilter2 <- function(peaklist, nsam = 2, p = 0.05){
  
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
  cat("done.");
  
  ##(2) prepare new peaklist
  cat("\n(2) Performing the first filtering...");
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
  
  ##(3) statistical test
  cat("\n(3) Performing statistical test...");
  ## generate an empty matrix 
  aov.out <- vector("list", dim(peaklist)[1])
  TukeyHSD.out <- vector("list", dim(peaklist)[1])
  ## get each variable names (original row numbers)
  var_peaks <- names(peaklist_new )[1:dim(peaklist)[1]] 
  ## perform test
  aov.out <- lapply(var_peaks, function(x) {
    lm(substitute(i ~ group, list(i = as.name(x))), data = peaklist_new)}) 
  ## extracting adjusted p-value
  TukeyHSD.out <- lapply(aov.out, function(x) TukeyHSD(aov(x))$group[, 4])
  output <- data.frame(matrix(unlist(TukeyHSD.out), ncol = 3, byrow = TRUE))
  colnames(output) <- c("CB", "BD", "CD")
  cat("done.");
  
  ##(4) perform the first filter
  cat("\n(4) Performing the first filtering...");
  c_index <- which(output$CB <= p & output$CD <= p)
  d_index <- which(output$BD <= p & output$CD <= p)
  peaklistB <- peaklist[peaklist$B > 0, ]
  peaklistC <- peaklist[c_index, ]
  peaklistD <- peaklist[d_index, ]
  exp.B <- cbind.data.frame(mz = peaklistB$mz, intensity = peaklistB[, (13 + each_exp[1])],
                            rt = peaklistB$rt)
  exp.C <- cbind.data.frame(mz = peaklistC$mz, intensity = peaklistC[, (13 + sum(each_exp[1:2]))],
                            rt = peaklistC$rt)
  exp.D <- cbind.data.frame(mz = peaklistD$mz, intensity = peaklistD[, (13 + sum(each_exp[1:3]))],
                            rt = peaklistD$rt)
  exp_list <- list(exp.B = exp.B, exp.C = exp.C, exp.D = exp.D)
  cat("done.");
  return(exp_list)
}

