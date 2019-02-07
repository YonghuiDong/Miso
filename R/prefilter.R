#' @title Prefilter
#' @description prefiltering isotopically labeled analytes according to the experiment design.
#' @param xset xcms object.
#' @param subgroup subset the xcms groups. The name should be the same as in phboData$class. default = NULL, which means no subset will be performed.
#' @param unlabel specify which is unlabeled group.
#' @param reps if there are replicates in the sample.
#' @param p p-value threshold, default value = 0.05
#' @param folds fold change threshold, default value = 10
#' @importFrom xcms peakTable
#' @importFrom stats lm
#' @importFrom stats aov
#' @importFrom stats TukeyHSD
#' @return a filtered peaklist
#' @export
#' @examples
#' data(lcms)
#' explist <- prefilter(lcms, subgroup = c("B", "C", "D"), unlabel = "B")

prefilter <- function(xset, subgroup = NULL, unlabel = NULL,  reps = TRUE, p = 0.05, folds = 10){

  peakTable <- NULL
  ##(1) check input
  cat("\n(1) Checking input parameters...");
  ## check object type
  if(class(xset) != "xcmsSet") {stop("the input object is not an xcmsSet object")}
  ## check xset phenoData
  pheno_levels <- levels(xset@phenoData$class)
  if(length(pheno_levels) < 2) {stop("at least two sample groups should be included")}
  ## check subsetgroup
  if(is.null(subgroup) == FALSE & all(subgroup %in% pheno_levels) == FALSE)
  {stop("selected subgroup(s) do not exist in your data")}
  ## check unlabel
  if(length(unlabel) == 0) {stop("unlabeled group does not exist")}
  if(unlabel %in% levels(xset@phenoData$class) == FALSE) {stop("unlabled group does not exist")}
  ## check reps
  if(reps %in% c(TRUE, FALSE) == FALSE) {stop("specify if replicates are contained in your sample: TURE/FALSE")}
  ## check p value
  if(!is.numeric(p)){stop("invalid calss of p-value threshold: not numeric")}
  if(p > 1 | p < 0){stop("invalid p-value, the value should between 0 and 1")}
  cat("done");

  ##(2) prepare new peaklist
  ## (2.1) subset peak
  ## get peaktable
  peak <- peakTable(xset)
  ## the first 7 rows are always the same in any aligned xcms peak table
  A = peak[, c(-1:-(7 + length(pheno_levels)))]
  B = t(A)
  C = cbind.data.frame(B, group = xset@phenoData$class)
  ## only select subgroups if subgroup is not NULL
  if(is.null(subgroup) == TRUE) {
    peaklist_new = C
  } else {
    peaklist_new <- C[(C$group %in% subgroup),]
    ## drop factors
    peaklist_new$group <- factor(peaklist_new$group)
  }

  ##(3) calculating fold change
  cat("\n(2) Calculating fold change...");
  fold <- fold(peaklist_new)
  ## add fold infor in peak
  peak <- cbind(fold, peak)
  cat("done");

  ##(4) statistical test
  cat("\n(3) Performing statistical test...");
  if (reps == TRUE) {
    stat <- getp(peaklist_new)
  } else {
    stat <- getp0(peaklist_new)
  }

  cat("done");

  ##(5) perform the first filter
  cat("\n(4) Performing the first filtering...");

  ## extract the index in which either p is less than the p-value threshold or fold change
  ## is higher than the fold change threshold
  ##(5.1) creat an ampty list that will store the filtered data
  n_index <- length(levels(peaklist_new$group))
  exp_list <- vector("list", length = n_index)
  names(exp_list) <- levels(peaklist_new$group)

  ##(5.2) data filtering
  for (i in 1 : n_index){
   group_name = levels(peaklist_new$group)[i]
   # subset data according to group
   data_p = as.data.frame(stat[, grepl(group_name, colnames(stat))])
   data_f = as.data.frame(fold[,grepl(paste(group_name, "_", sep = ""), colnames(fold))])
   index_p <-  which(apply(data_p < p, 1, all) == TRUE) # get p index
   index_f <- which(apply(data_f >= folds, 1, all) == TRUE) # get fold index
   index <- sort(unique(c(index_p, index_f)))
   ## first filter with p value and fold change
   peak_filter1 <- peak[index, ]
   ## second filter with detected peak number, and select specific columns
   exp_list[[i]] <- peak_filter1[peak_filter1[group_name] > 0,][, c("mz", "rt", paste("mean_", group_name, sep = ""))]
   colnames(exp_list[[i]]) <- c("mz", "rt", "intensity")
  }

  ##(5.3) refilter the unlabeled group
  exp_list[[unlabel]] <- peak[peak[unlabel] > 0, ][, c("mz", "rt", paste("mean_", unlabel, sep = ""))]
  colnames(exp_list[[unlabel]]) <- c("mz", "rt", "intensity")
  exp_list$stat <- stat
  exp_list$fold <- fold
  cat("done");
  return(exp_list)
}

