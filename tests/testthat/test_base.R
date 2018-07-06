context("Output check")

test_that("varImp functions", {
  library(party)
  library(measures)
  
  iris.cf <- cforest(Species ~ ., data = iris,control = cforest_unbiased(mtry = 2, ntree = 50))
  varImpAUC(object = iris.cf, method = "ovo")
  # Erweitern auf beliebige Maße?
  
  # binary case
  readingSkills.cf <- cforest(score ~ ., data = readingSkills, 
    control = cforest_unbiased(mtry = 2, ntree = 50))
  varImpAUC(object = readingSkills.cf, method = "ovo")
  
  iris2 = iris
  iris2$Species = factor(iris$Species == "versicolor")
  iris.cf <- cforest(Species ~ ., data = iris2,control = cforest_unbiased(mtry = 2, ntree = 50))
  set.seed(123)
  a = varImpAUC(object = iris.cf, method = "ova")
  set.seed(123)
  b = varImpAUC(object = iris.cf, method = "ovo")
  set.seed(123)
  # c = varimpAUC(object = iris.cf) # current party implementation is wrong!
  expect_equal(a,b)
  # expect_equal(b,c)
  
  # multiclass case
  iris.cf <- cforest(Species ~ ., data = iris, control = cforest_unbiased(mtry = 2, ntree = 50))
  set.seed(123)
  a = varImpAUC(object = iris.cf, method = "ovo")
  set.seed(123)
  b = varImpAUC(object = iris.cf, method = "ova")
  set.seed(123)
  c = varImpACC(object = iris.cf)
  
  expect_true(all(!is.na(a)))
  expect_true(all(!is.na(b)))
  expect_true(all(!is.na(c)))
})
