library(ranger)
library(measures)

rg.iris <- ranger(Species ~ ., data = iris, keep.inbag = TRUE, probability = TRUE)
model = rg.iris
data = iris
target = "Species"

varImp = function(model, data, target, nperm = 1, measure = "multiclass.Brier") {
  pred_cols = which(colnames(data) != target)
  num.trees = model$num.trees
  inbag = do.call(cbind, model$inbag.counts)
  pred_levels = levels(data[, target])
  truth = data[, target]
  
  # Calculate original performance
  old_predis = predict(model, data = data, predict.all = TRUE)$predictions
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

      # old slow version
      # preds = predict(model, data = data_new, predict.all = TRUE)
      # predis = preds$predictions[, , i]
      
      # fast version
      model2 <- swap_trees(model, 1, i)
      predis <- predict(model2, data = data_new, num.trees = 1)$predictions

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

# Old version
# 0.07379483 0.01249709 0.58453166 0.55452196
# User      System verstrichen 
# 38.940       7.796      31.104 

# New version
# 0.07558051 0.01276306 0.58440709 0.55576188
# User      System verstrichen 
# 10.840       1.676      11.173 

library(profvis)
profvis({vimp_brier <- varImp(model, data, target)})
vimp_brier

