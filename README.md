# varImp
Random forest variable importance for arbitrary measures of the [measures](https://github.com/mlr-org/measures) package, which contains the biggest collection of measures for regression and classification in R. 

## Installation
The development version

    devtools::install_github("mlr-org/measures")
    devtools::install_github("PhilippPro/varImp")
    iris.cf <- cforest(Species ~ ., data = iris, control = cforest_unbiased(mtry = 2, ntree = 50))
    varImp(object = iris.cf, measure = "multiclass.Brier")
    varImpACC(object = iris.cf)
    # Two classes:
    iris2 = iris
    iris2$Species = factor(iris$Species == "versicolor")
    iris.cf = cforest(Species ~ ., data = iris2,control = cforest_unbiased(mtry = 2, ntree = 50))
    varImpAUC(object = iris.cf)