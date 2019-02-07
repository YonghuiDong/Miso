#' @title get artifical p-values
#' @description get artifical p-values for groups without replicates
#' @param x dataframe
#' @return a data frame
#' @export
#' @examples
#' dat = as.data.frame(matrix(runif(2*3), ncol = 2, nrow = 3))
#' dat$group = rep(LETTERS[2:4], 1)
#' out <- getp0(dat)

getp0 <- function(x){
  f <- as.factor(x$group)
  ## rm group info in order to calculate colMeans 
  j <- combn(levels(f), 2)
  output <- matrix(0, nrow = dim(x)[2] - 1, ncol = dim(j)[2])
  colnames(output) <- paste(j[1,], "-", j[2,], sep = '')
  output <- as.data.frame(output)
  return(output)
}
