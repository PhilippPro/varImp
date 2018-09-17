library(ranger)
library(measures)

rg.iris <- ranger(Species ~ ., data = iris, keep.inbag = TRUE, probability = TRUE)
object = rg.iris
data = iris
target = "Species"

varImp = function(object, data, target, nperm = 1, measure = "multiclass.Brier") {
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

  # todo: Put the predictions internal/only for one tree, to make them much faster
 
  res_new = matrix(NA, num.trees, length(pred_cols))
  
  for(j in pred_cols) {
    for(i in 1:num.trees) { 
      print(paste("column", j, "tree", i))
      data_new = data[inbag[,i] == 0,]
      data_new[,j] = sample(data_new[,j], replace = FALSE)
      preds = predict(object, data = data_new, predict.all = TRUE)
      predis = preds$predictions[, , i]
      colnames(predis) = pred_levels
      truth = data_new[, target]
      perf_new = do.call(measure, list(predis, truth))
      # minSign = ifelse(MEASUREMINIMIZE, 1, -1)
      res_new[i, j] = perf_new 
    }
  }
  res = res_new - res_old
  
  return(colMeans(res))
}

library(profvis)
profvis({vimp_brier <- varImp(object, data, target)})
vimp_brier

