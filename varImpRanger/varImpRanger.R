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
  # Evtl. viel sinnvoller/schneller die Predictions am Ende zu aggregieren und dann das measure zu berechnen; hier muss man überlegen
 
  erg = matrix(NA, num.trees, length(pred_cols))
  
  for(j in pred_cols) {
    for(i in 1:num.trees) { 
      print(paste("column", j, "tree", i))
      data_new = data[inbag[,i] == 0,]
      preds = predict(object, data = data_new, predict.all = TRUE)
      predis = preds$predictions[, , 1]
      colnames(predis) = pred_levels
      truth = data_new[, target]
      perf_old = do.call(measure, list(predis, truth))
      
      data_new[,j] = sample(data_new[,j], replace = FALSE)
      preds = predict(object, data = data_new, predict.all = TRUE)
      predis = preds$predictions[, , 1]
      colnames(predis) = pred_levels
      truth = data_new[, target]
      perf_new = do.call(measure, list(predis, truth))
      # minSign = ifelse(MEASUREMINIMIZE, 1, -1)
      erg[i, j] = perf_new - perf_old
    }
  }
  return(colMeans(erg))
}

library(profvis)
profvis({vimp_brier <- varImp(object, data, target)})
vimp_brier

