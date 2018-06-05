# varImp
Variable Importance for arbitrary measures

## Installation
The development version

    devtools::install_github("mlr-org/measures")
    devtools::install_github("PhilippPro/varImp")
    iris.cf <- cforest(Species ~ ., data = iris, control = cforest_unbiased(mtry = 2, ntree = 50))
    varImp(object = iris.cf, measure = "multiclass.Brier")