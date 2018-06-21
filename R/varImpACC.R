#' varImpACC
#' 
#' Computes the variable importance regarding the accuracy (ACC). 
#'
#' @param object an object as returned by cforest.
#' @param mincriterion the value of the test statistic or 1 - p-value that must be exceeded in order to include a 
#' split in the computation of the importance. The default mincriterion = 0 guarantees that all splits are included.
#' @param conditional the value of the test statistic or 1 - p-value that must be exceeded in order to include a split 
#' in the computation of the importance. The default mincriterion = 0 guarantees that all splits are included.
#' @param threshold the threshold value for (1 - p-value) of the association between the variable of interest and a 
#' covariate, which must be exceeded inorder to include the covariate in the conditioning scheme for the variable of 
#' interest (only relevant if conditional = TRUE). A threshold value of zero includes all covariates.
#' @param nperm the number of permutations performed.
#' @param OOB a logical determining whether the importance is computed from the out-of-bag sample or the learning 
#' sample (not suggested).
#' @param pre1.0_0 Prior to party version 1.0-0, the actual data values were permuted according to the original 
#' permutation importance suggested by Breiman (2001). Now the assignments to child nodes of splits in the variable 
#' of interest are permuted as described by Hapfelmeier et al. (2012), which allows for missing values in the 
#' explanatory variables and is more efficient wrt memory consumption and computing time. This method does not 
#' apply to conditional variable importances.
#'
#' @return vector with computed permutation importance for each variable
#' @export
#'
#' @examples
#' data(iris)
#' iris2 = iris
#' iris2$Species = factor(iris$Species == "versicolor")
#' iris.cf <- cforest(Species ~ ., data = iris2,control = cforest_unbiased(mtry = 2, ntree = 50))
#' set.seed(123)
#' a = varImpACC(object = iris.cf)
#' 
varImpACC <- function (object, mincriterion = 0, conditional = FALSE, threshold = 0.2, 
  nperm = 1, OOB = TRUE, pre1.0_0 = conditional) { 
  response <- object@responses
  input <- object@data@get("input")
  xnames <- colnames(input)
  inp <- initVariableFrame(input, trafo = NULL)
  y <- object@responses@variables[[1]]
  if (length(response@variables) != 1) 
    stop("cannot compute variable importance measure for multivariate response")
  if (conditional || pre1.0_0) {
    if (!all(complete.cases(inp@variables))) 
      stop("cannot compute variable importance measure with missing values")
  }
  CLASS <- all(response@is_nominal)
  ORDERED <- all(response@is_ordinal)
  
  error <- function(x, oob) {
    xoob <- t(sapply(x, function(x) x))[oob,]
    yoob <- y[oob]
    return(measures::ACC(xoob, yoob))
  }
  w <- object@initweights
  if (max(abs(w - 1)) > sqrt(.Machine$double.eps)) 
    warning(sQuote("varimp"), " with non-unity weights might give misleading results")
  perror <- matrix(0, nrow = nperm * length(object@ensemble), ncol = length(xnames))
  colnames(perror) <- xnames
  for (b in 1:length(object@ensemble)) {
    tree <- object@ensemble[[b]]
    if (OOB) {
      oob <- object@weights[[b]] == 0
    } else {
      oob <- rep(TRUE, length(xnames))
    }
    p <- party_intern(tree, inp, mincriterion, -1L, fun = "R_predict") 
    eoob <- error(p, oob)
    for (j in unique(varIDs(tree))) {
      for (per in 1:nperm) {
        if (conditional || pre1.0_0) {
          tmp <- inp
          ccl <- create_cond_list(conditional, threshold, 
            xnames[j], input)
          if (is.null(ccl)) {
            perm <- sample(which(oob))
          }
          else {
            perm <- conditional_perm(ccl, xnames, input, 
              tree, oob)
          }
          tmp@variables[[j]][which(oob)] <- tmp@variables[[j]][perm]
          p <- party_intern(tree, tmp, mincriterion, -1L, fun = "R_predict") 
        } else {
          p <- party_intern(tree, inp, mincriterion, as.integer(j), fun = "R_predict") 
        }
        perror[(per + (b - 1) * nperm), j] <- - (error(p,oob) - eoob)
      }
    }
  }
  perror <- as.data.frame(perror)
  return(MeanDecrease = colMeans(perror, na.rm = TRUE))
}

