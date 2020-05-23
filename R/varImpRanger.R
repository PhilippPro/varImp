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
#' \dontrun{
#' library(ranger)
#' iris.rg = ranger(Species ~ ., data = iris, keep.inbag = TRUE, probability = TRUE)
#' vimp.ranger = varImpRanger(object = iris.rg, data = iris, target = "Species")
#' vimp.ranger
#' }
varImpRanger = function(object, data, target, nperm = 1, measure = "multiclass.Brier") {
  # Some tests
  if(!("ranger" %in% class(object)))
    stop("Object is not a 'ranger' model")
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
  for(i in 1:num.trees)
    res_old[i] = do.call(measure, list(old_predis[inbag[,i] == 0, , i],  truth[inbag[,i] == 0]))
  
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
  
  return(colMeans(res))
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