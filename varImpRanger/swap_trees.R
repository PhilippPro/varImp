library(ranger)

swap_trees <- function(rf, i, j) {
  res <- rf
  
  res$forest$child.nodeIDs[[i]] <- rf$forest$child.nodeIDs[[j]]
  res$forest$split.varIDs[[i]] <- rf$forest$split.varIDs[[j]]
  res$forest$split.values[[i]] <- rf$forest$split.values[[j]]
  
  res$forest$child.nodeIDs[[j]] <- rf$forest$child.nodeIDs[[i]]
  res$forest$split.varIDs[[j]] <- rf$forest$split.varIDs[[i]]
  res$forest$split.values[[j]] <- rf$forest$split.values[[i]]
  
  if (!is.null(rf$forest$terminal.class.counts)) {
    res$forest$terminal.class.counts[[i]] <- rf$forest$terminal.class.counts[[j]]
    res$forest$terminal.class.counts[[j]] <- rf$forest$terminal.class.counts[[i]]
  }
  
  res
}

rf <- ranger(Species ~., iris, num.trees = 5, probability = TRUE)
rf2 <- swap_trees(rf, 1, 4)
pred_tree4 <- predict(rf2, iris, num.trees = 1)$predictions
# Check
pred_all <- predict(rf, iris, predict.all = TRUE)$predictions
pred_tree4 == pred_all[,,4]