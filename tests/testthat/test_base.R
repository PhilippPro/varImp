context("Output check")

test_that("varImp functions", {
  library(party)
  library(measures)
  
  # regression
  readingSkills.cf = cforest(score ~ ., data = readingSkills, 
    control = cforest_unbiased(mtry = 2, ntree = 50))
  varImp(object = readingSkills.cf, measure = "MSE")
  
  # Erweitern auf beliebige Maße?
  
  # binary case
  iris2 = iris
  iris2$Species = factor(iris$Species == "versicolor")
  iris.cf = cforest(Species ~ ., data = iris2,control = cforest_unbiased(mtry = 2, ntree = 50))
  set.seed(123)
  a = varImpAUC(object = iris.cf)
  expect_true(all(!is.na(a)))

  # expect_equal(b,c)
  d = varImp(object = iris.cf, measure = "Brier", positive = "FALSE")
  e = varImp(object = iris.cf, measure = "ACC")
  f = varImp(object = iris.cf, measure = "MMCE")
  expect_true(all(!is.na(d)))
  expect_true(all(!is.na(e)))
  expect_true(all(!is.na(f)))
  
  # multiclass case
  iris.cf = cforest(Species ~ ., data = iris, control = cforest_unbiased(mtry = 2, ntree = 50))
  set.seed(123)
  a = varImp(object = iris.cf, measure = "multiclass.AU1P")
  set.seed(123)
  b = varImpACC(object = iris.cf)
  
  expect_true(all(!is.na(a)))
  expect_true(all(!is.na(b)))

  # ranger
  library(ranger)
  iris.rg = ranger(Species ~ ., data = iris, keep.inbag = TRUE, probability = TRUE)
  vimp.ranger = varImpRanger(object = iris.rg, data = iris, target = "Species")
  expect_true(is.numeric(vimp.ranger))
  expect_true(all(!is.na(vimp.ranger)))
})
