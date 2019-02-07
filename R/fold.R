#' @title fold change
#' @description calculate fold change among different samples.
#' @param x dataframe
#' @importFrom utils combn
#' @return a dataframe with mean values and fold changes
#' @export
#' @examples
#' vars <- 1000
#' samples <- 50
#' groups <- 3
#' dat <- replicate(vars, runif(n = samples))
#' f <- rep_len(1:groups, samples)
#' f <- LETTERS[f]
#' dat <- data.frame(dat, group = f)
#' ret <- fold(dat)

 fold <- function(x){
  f <- x$group
  i <- split(1:nrow(x), f)
  ## rm group info in order to calculate colMeans 
  x$group = NULL
  mean_int <- sapply(i, function(i){colMeans(x[i, ])})
  x <- t(mean_int)
  j <- combn(levels(f), 2)
  f_change1 <- x[j[1,],] / x[j[2,],]
  f_change2 <- x[j[2,],] / x[j[1,],]
  ## remove NaN in f_change Matrix
  f_change <- rbind(f_change1, f_change2)
  f_change[is.nan(f_change)] <- 0
  rownames(f_change) <- c(paste(j[1,], "_", j[2,], sep = ''), paste(j[2,], "_", j[1,], sep = '') )
  f_change <- as.data.frame(t(f_change))
  colnames(mean_int) <- paste("mean_", colnames(mean_int), sep = "")
  ret <- cbind.data.frame(mean_int, f_change)
  ret
}
