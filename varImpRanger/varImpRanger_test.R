#' varImpRanger
#' 
#' Computes the variable importance for ranger models and for arbitrary measures from the 'measures' package.
#'
#' @param object An object as returned by cforest. \code{\link[ranger]{ranger}} with option \code{keep.inbag = TRUE}.
#' @param data Original data that was used for training the random forest.
#' @param target Target variable as used in the trained model.
#' @param nperm The number of permutations performed.
#' @param measure The name of the measure of the 'measures' package that should be used for the variable importance calculation.
#'
#'
#' @importFrom stats predict
#' @return Vector with computed permutation importance for each variable.
#' @export
#'
#' @examples
#' library(ranger)
#' iris.rg = ranger(Species ~ ., data = iris, keep.inbag = TRUE, probability = TRUE)
#' vimp.ranger = varImpRanger(object = iris.rg, data = iris, target = "Species")
#' vimp.ranger
varImpRanger_test = function(object, data, target, nperm = 1, measure = "multiclass.Brier") {
  # Some tests
  if(!("ranger" %in% class(object)))
    stop("Object is not a 'ranger' model")
  if(!(object$treetype %in% c("Regression", "Probability estimation")))
    stop("Object is not a regression or probability forest")
  measureList = listAllMeasures()
  if (!(measure %in% measureList[, 1]))
    stop("measure should be a measure of the measures package")
  measure.minimize = measureList$minimize[measureList[,1] == measure]
  
  pred_cols = which(colnames(data) != target)
  num.trees = object$num.trees
  inbag = do.call(cbind, object$inbag.counts)
  pred_levels = levels(data[, target])
  truth = data[, target]
  
  # Calculate original performance
  old_predis = predict(object, data = data, predict.all = TRUE)$predictions
  colnames(old_predis) = pred_levels
  res_old = numeric(num.trees)
  for(i in 1:num.trees){
    if(object$treetype == "Regression")
      res_old[i] = do.call(measure, list(old_predis[inbag[,i] == 0,  i],  truth[inbag[,i] == 0]))
    if(object$treetype == "Probability estimation")
      res_old[i] = do.call(measure, list(old_predis[inbag[,i] == 0, , i],  truth[inbag[,i] == 0]))
  }
  
  # Calculate permuted performance
  res_new = matrix(NA, num.trees, length(pred_cols))
  for(j in pred_cols) {
    for(i in 1:num.trees) { 
      print(paste("column", j, "tree", i))
      data_new = data[inbag[,i] == 0,]
      data_new[,j] = sample(data_new[,j], replace = FALSE)
      object2 = swap_trees(object, 1, i)
      predis = predict(object2, data = data_new, num.trees = 1)$predictions
      
      colnames(predis) = pred_levels
      truth = data_new[, target]
      perf_new = do.call(measure, list(predis, truth))
      
      res_new[i, j] = perf_new 
    }
  }
  minSign = ifelse(measure.minimize, 1, -1)
  res = minSign * (res_new - res_old)
  
  return(list(res = res, vimp = colMeans(res)))
}
swap_trees = function(rf, i, j) {
  res = rf
  
  res$forest$child.nodeIDs[[i]] = rf$forest$child.nodeIDs[[j]]
  res$forest$split.varIDs[[i]] = rf$forest$split.varIDs[[j]]
  res$forest$split.values[[i]] = rf$forest$split.values[[j]]
  
  res$forest$child.nodeIDs[[j]] = rf$forest$child.nodeIDs[[i]]
  res$forest$split.varIDs[[j]] = rf$forest$split.varIDs[[i]]
  res$forest$split.values[[j]] = rf$forest$split.values[[i]]
  
  if (!is.null(rf$forest$terminal.class.counts)) {
    res$forest$terminal.class.counts[[i]] = rf$forest$terminal.class.counts[[j]]
    res$forest$terminal.class.counts[[j]] = rf$forest$terminal.class.counts[[i]]
  }
  res
}

devtools::install_github("mlr-org/measures")
devtools::install_github("PhilippPro/varImp")
library(ranger)
library(varImp)

iris.rg = ranger(Species ~ ., data = iris, keep.inbag = TRUE, probability = TRUE)
vimp.ranger = varImpRanger_test(object = iris.rg, data = iris, target = "Species", measure = "multiclass.Brier") # MMCE geht nicht!

# binomial test
bin.test_greater = function(y) {
  x = sum(y > 0)
  n = length(y)
  binom.test(x, n, p = 0.5, alternative = "greater")$p.value
}

bin.test_two.sided = function(y) {
  x = sum(y > 0)
  n = length(y)
  binom.test(x, n, p = 0.5, alternative = "two.sided")$p.value
}

# z-Scores
p_value = function(x) 1 - pnorm(mean(x) / (sd(x)/sqrt(length(x))))

vimp.ranger$vimp
apply(vimp.ranger$res < 0, 2, mean)
round(apply(vimp.ranger$res, 2, bin.test), 3)
apply(vimp.ranger$res, 2, p_value)

# Simulation

sim_data = function(n = 20) {
  x1 = runif(n,5,95)
  x2 = runif(n,5,95)
  x3 = rbinom(n,1,.5)
  
  b0 = 17
  b1 = 0
  b2 = 0.037
  b3 = -5.2
  sigma = 1.4
  
  eps = rnorm(n, 0, sigma)
  y = b0 + b1*x1 + b2*x2  + b3*x3 + eps
  X = cbind(x1, x2, x3)
  return(data.frame(X, y))
}

sim = sim_data(n = 1000)
mod = ranger(y ~ ., data = sim, num.trees = 1000, keep.inbag = TRUE)
vimp.ranger = varImpRanger_test(object = mod, data = sim, target = "y", measure = "MSE")
apply(vimp.ranger$res < 0, 2, mean)
# n = 100, 1000 Baeume: 0.429 0.243 0.000
# n = 1000, 1000 Baeume: 0.411 0.006 0.000
# n = 10000, 1000 Baeume: 0.496 0.003 0.000
# n = 100000, 1000 Baeume: 0.469 0.002 0.000

round(apply(vimp.ranger$res, 2, bin.test_greater), 3)
# n = 100, 1000 Baeume: 0.98 0.00 0.00
# n = 1000, 1000 Baeume: 0.747 0.000 0.000
# n = 10000, 1000 Baeume: 0.803 0.000 0.000
# Scheint dann öfter > 0 zu sein, wenn nicht relevant! (warum?)

round(apply(vimp.ranger$res, 2, bin.test_two.sided), 3)
# n = 100, 1000 Baeume: 0.046 0.000 0.000
# n = 1000, 1000 Baeume: 0.548 0.000 0.000
# n = 10000, 1000 Baeume: 0.018 0.000 0.000

apply(vimp.ranger$res, 2, p_value)
# n = 100, 1000 Baeume: 1.895139e-05 0.000000e+00 0.000000e+00
# n = 1000, 1000 Baeume: 0.09871624 0.00000000 0.00000000
# n = 10000, 1000 Baeume: 0.695337 0.000000 0.000000

summary(lm(y~., data = sim))
