#' @title get p-values
#' @description get p-values from Post Hoc analysis
#' @param dat peaklist
#' @importFrom stats as.formula
#' @importFrom stats formula
#' @importFrom stats terms
#' @importFrom stats terms.formula
#' @return a data frame
#' @export
#' @examples
#' dat = as.data.frame(matrix(runif(2*300), ncol = 2, nrow = 300))
#' dat$group = rep(LETTERS[2:4], 100)
#' out <- getp(dat)

getp <- function(dat){
  ##(1) format the data
  dat$group <- factor(dat$group , levels = c("D", "C", "B"))
  colnames(dat)[-ncol(dat)] <- paste("V", names(dat)[-ncol(dat)], sep = "")
  ##(2) aov
  response_names <- names(dat)[-ncol(dat)]
  form <- as.formula(sprintf("cbind(%s) ~ group", toString(response_names)))
  fit <- do.call("aov", list(formula = form, data = quote(dat)))
  aov_hack <- fit
  aov_hack[c("coefficients", "fitted.values")] <- NULL
  aov_hack[c("contrasts", "xlevels")] <- NULL
  attr(aov_hack$model, "terms") <- NULL
  class(aov_hack) <- c("aov", "lm")

  ##(3) post hoc
  N <- length(response_names)
  result <- vector("list", N)
  for (i in 1:N) {
    ## change response variable in the formula
    aov_hack$call[[2]][[2]] <- as.name(response_names[i])
    ## change residuals
    aov_hack$residuals <- fit$residuals[, i]
    ## change effects
    aov_hack$effects <- fit$effects[, i]
    ## change "terms" object and attribute
    old_tm <- terms(fit)  ## old "terms" object in the model
    old_tm[[2]] <- as.name(response_names[i])  ## change response name in terms
    new_tm <- terms.formula(formula(old_tm))  ## new "terms" object
    aov_hack$terms <- new_tm  ## replace `aov_hack$terms`
    ## replace data in the model frame
    aov_hack$model[1] <- data.frame(fit$model[[1]][, i])
    names(aov_hack$model)[1] <- response_names[i]
    ## run `TukeyHSD` on `aov_hack`
    result[[i]] <- TukeyHSD(aov_hack)$group[, 4]
  }
  ##(4) output
  output <- data.frame(matrix(unlist(result), ncol = 3, byrow = TRUE))
  colnames(output) <- c("p_CD", "p_BD", "p_CB")
  output
}


