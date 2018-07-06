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
  return(varImp(object, mincriterion = mincriterion, conditional = conditional, threshold = threshold, nperm = nperm, 
    OOB = OOB, pre1.0_0 = pre1.0_0, measure = "ACC"))
}

