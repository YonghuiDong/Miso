#' @title fold change
#' @description calculate fold change among different samples.
#' @param x matrix data, type: matrix
#' @param f groups
#' @importFrom utils combn
#' @return a dataframe with fold changes
#' @export
#' @examples
#' vars <- 1000
#' samples <- 50
#' groups <- 3
#' dat <- replicate(vars, runif(n = samples))
#' f <- rep_len(1:groups, samples) + 1
#' f <- LETTERS[f]
#' ret <- fold(dat, f)

 fold <- function(x, f){
  f <- factor(f, levels = c("D", "C", "B"))
  i <- split(1:nrow(x), f)
  mean_int <- sapply(i, function(i){colMeans(x[i,])})
  x <- t(mean_int)
  j <- combn(levels(f), 2)
  f_change <- x[j[1,],] / x[j[2,],]
  ## remove NaN in f_change Matrix
  f_change[is.nan(f_change)] <- 0
  rownames(f_change) <- paste('F_', j[1,], j[2,], sep = '')
  f_change <- as.data.frame(t(f_change))
  ret <- cbind.data.frame(mean_int, f_change, F_CD = 1 / f_change$F_DC)
  ret
}
