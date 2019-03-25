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
bin.test_greater = function(y, zeros = FALSE) {
  #if(!(zeros))
  #  y = y[y!=0] #  Diese Loesung ist nicht wirklich sauber
  x = sum(y > 0)
  n = length(y)
  binom.test(x, n, p = 0.5, alternative = "greater")$p.value
}

bin.test_two.sided = function(y, zeros = FALSE) {
  #if(!(zeros))
  #  y = y[y!=0] #  Diese Loesung ist nicht wirklich sauber
  x = sum(y > 0)
  n = length(y)
  binom.test(x, n, p = 0.5, alternative = "two.sided")$p.value
}

# wilcoxon test
wilcoxon.test = function(y, zeros = FALSE) {
  #if(!(zeros))
  #  y = y[y!=0] #  Diese Loesung ist nicht wirklich sauber
  wilcox.test(y, mu = 10^-10, alternative = "greater")$p.value
}

# z-Scores
p_value = function(y, zeros = FALSE) {
  #if(!(zeros))
  #  y = y[y!=0] #  Diese Loesung ist nicht wirklich sauber
  t.test(y, alternative = "greater", mu = 0)$p.value #  1 - pnorm(mean(x) / (sd(x)/sqrt(length(x))))
} 

vimp.ranger$vimp
apply(vimp.ranger$res < 0, 2, mean)
round(apply(vimp.ranger$res, 2, bin.test_greater), 3)
round(apply(vimp.ranger$res, 2, bin.test_two.sided), 3)
round(apply(vimp.ranger$res, 2, wilcoxon.test), 3)
round(apply(vimp.ranger$res, 2, p_value), 3)

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
mod = ranger(y ~ ., data = sim, num.trees = 10000, keep.inbag = TRUE)
vimp.ranger = varImpRanger_test(object = mod, data = sim, target = "y", measure = "MSE")
apply(vimp.ranger$res[1:10000,] < 0, 2, mean)
# n = 100, 1000 Baeume: 0.429 0.243 0.000
# n = 1000, 1000 Baeume: 0.411 0.006 0.000
# n = 10000, 1000 Baeume: 0.496 0.003 0.000
# n = 100000, 1000 Baeume: 0.469 0.002 0.000

round(apply(vimp.ranger$res, 2, bin.test_greater), 3)
# n = 100, 1000 Baeume: 0.98 0.00 0.00
# n = 1000, 1000 Baeume: 0.747 0.000 0.000
# n = 10000, 1000 Baeume: 0.803 0.000 0.000
# Scheint dann öfter > 0 zu sein, wenn nicht relevant? Ist dem überhaupt so? -> öfters simulieren...

round(apply(vimp.ranger$res, 2, bin.test_two.sided), 3)
# n = 100, 1000 Baeume: 0.046 0.000 0.000
# n = 1000, 1000 Baeume: 0.548 0.000 0.000
# n = 10000, 1000 Baeume: 0.018 0.000 0.000

round(apply(vimp.ranger$res, 2, wilcoxon.test), 3)

round(apply(vimp.ranger$res, 2, p_value), 3)
# n = 100, 1000 Baeume: 1.895139e-05 0.000000e+00 0.000000e+00
# n = 1000, 1000 Baeume: 0.09871624 0.00000000 0.00000000
# n = 10000, 1000 Baeume: 0.695337 0.000000 0.000000

a = summary(lm(y~., data = sim))

# Distribution of the p_values under the Null-Hypothesis

n_exp = 1000
p_matrix = array(NA, dim = c(n_exp, 7, 3))
resis = list()

for(i in 1:100) {
  print(i)
  set.seed(i)
  sim = sim_data(n = 100)
  mod = ranger(y ~ ., data = sim, num.trees = 10000, keep.inbag = TRUE)
  invisible(capture.output(vimp.ranger <- varImpRanger_test(object = mod, data = sim, target = "y", measure = "MSE")))
  resis[[i]] = vimp.ranger$res
  p_matrix[i, 1, ] = apply(resis[[i]] < 0, 2, mean)
  p_matrix[i, 2, ] = apply(resis[[i]] > 0, 2, mean)
  p_matrix[i, 3, ] = round(apply(resis[[i]], 2, bin.test_greater), 10)
  p_matrix[i, 4, ] = round(apply(resis[[i]], 2, bin.test_two.sided), 10)
  p_matrix[i, 5, ] = round(apply(resis[[i]], 2, wilcoxon.test), 10)
  p_matrix[i, 6, ] = round(apply(resis[[i]], 2, p_value), 10)
  p_matrix[i, 7, ] = summary(lm(y~., data = sim))$coefficients[2:4,4]

  save(resis, p_matrix, file = "./varImpRanger/simulation.RData")
}

load("./varImpRanger/simulation.RData")

par(mfrow = c(4,8))
mittel = numeric(100)
for(i in 1:100) {
  hist(resis[[i]][,1], main = i)
  print(mittel[i] <- mean(resis[[i]][,1]))
}

mean(mittel)
mean(mittel >0)
which(mittel > 0)
which(p_matrix[, 3, 1] <0.05)

sum(p_matrix[, 3, 1] >0.05, na.rm = T)
sum(p_matrix[, 4, 1] >0.05, na.rm = T)
sum(p_matrix[, 5, 1] >0.05, na.rm = T)
sum(p_matrix[, 6, 1] >0.05, na.rm = T)
sum(p_matrix[, 7, 1] >0.05, na.rm = T)

sum(resis[[i]][,1] > 0)
sum(resis[[i]][,1] <= 0)

load("./varImpRanger/simulation.RData")

# some Graphical analysis
hist(vimp.ranger$res[,1]) # zu viele exakt Null
a = vimp.ranger$res[,1]
sum(a==0)
mean(a)
mean(a[a!=0])
qqnorm(a)
qqline(a)
hist(a)
mean(a>0)
round(bin.test_greater(a), 3)
round(bin.test_greater(a, zeros = TRUE), 3)
round(p_value(a), 3)
round(p_value(a, zeros = TRUE), 3)

p_matrix[1:7, , ]

# meistens nur klare Trennungen bei  den anderen Tests...

for(i in 1:ncol(p_matrix))
  print(hist(p_matrix[,i, 1], main = i))




# Second simulation

n_exp = 1000
p_matrix = array(NA, dim = c(n_exp, 7, 3))
resis = list()

for(i in 1:100) {
  print(i)
  set.seed(i)
  sim = sim_data(n = 500)
  mod = ranger(y ~ ., data = sim, num.trees = 10000, keep.inbag = TRUE)
  invisible(capture.output(vimp.ranger <- varImpRanger_test(object = mod, data = sim, target = "y", measure = "MSE")))
  resis[[i]] = vimp.ranger$res
  p_matrix[i, 1, ] = apply(resis[[i]] < 0, 2, mean)
  p_matrix[i, 2, ] = apply(resis[[i]] > 0, 2, mean)
  p_matrix[i, 3, ] = round(apply(resis[[i]], 2, bin.test_greater), 10)
  p_matrix[i, 4, ] = round(apply(resis[[i]], 2, bin.test_two.sided), 10)
  p_matrix[i, 5, ] = round(apply(resis[[i]], 2, wilcoxon.test), 10)
  p_matrix[i, 6, ] = round(apply(resis[[i]], 2, p_value), 10)
  p_matrix[i, 7, ] = summary(lm(y~., data = sim))$coefficients[2:4,4]
  
  save(resis, p_matrix, file = "./varImpRanger/simulation2.RData")
}

load("./varImpRanger/simulation2.RData")

par(mfrow = c(4,8))
mittel = numeric(100)
for(i in 1:100) {
  hist(resis[[i]][,1], main = i)
  print(mittel[i] <- mean(resis[[i]][,1]))
}

mean(mittel)
mean(mittel >0)
which(mittel > 0)
which(p_matrix[, 3, 1] <0.05)

sum(p_matrix[, 3, 1] >0.05, na.rm = T)
sum(p_matrix[, 4, 1] >0.05, na.rm = T)
sum(p_matrix[, 5, 1] >0.05, na.rm = T)
sum(p_matrix[, 6, 1] >0.05, na.rm = T)
sum(p_matrix[, 7, 1] >0.05, na.rm = T)

sum(resis[[i]][,1] > 0)
sum(resis[[i]][,1] <= 0)

load("./varImpRanger/simulation.RData")

# some Graphical analysis
hist(vimp.ranger$res[,1]) # zu viele exakt Null
a = vimp.ranger$res[,1]
sum(a==0)
mean(a)
mean(a[a!=0])
qqnorm(a)
qqline(a)
hist(a)
mean(a>0)
round(bin.test_greater(a), 3)
round(bin.test_greater(a, zeros = TRUE), 3)
round(p_value(a), 3)
round(p_value(a, zeros = TRUE), 3)

p_matrix[1:7, , ]

# meistens nur klare Trennungen bei den anderen Tests...

for(i in 1:ncol(p_matrix))
  print(hist(p_matrix[,i, ]))

# 4 kann man wegschmeissen
# 3,5,6 könne man skalieren, so dass sie Varianz 1 haben?



# wilcoxon, binomial und silkes test
