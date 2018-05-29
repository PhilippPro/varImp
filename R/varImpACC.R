## mincriterion = 0 so that complete tree is evaluated; 
## regulate size of considered tree here via, e.g., mincriterion = 0.95
## or when building the forest in the first place via cforest_control(mincriterion = 0.95)

#' varImpACC
#'
#' @param object 
#' @param mincriterion 
#' @param conditional 
#' @param threshold 
#' @param nperm 
#' @param OOB 
#' @param pre1.0_0 
#'
#' @return vector with computed permutation importance for each variable
#' @export
#'
#' @examples
#' data(iris)
#' iris2 = iris
#' iris2$Species = factor(iris$Species == "versicolor")
#' iris.cf <- cforest(Species ~ ., data = iris2,control = cforest_unbiased(mtry = 2, ntree = 50))
#' set.seed(123)
#' a = varImpACC(object = iris.cf)
#' 
varImpACC <- function (object, mincriterion = 0, conditional = FALSE, 
  threshold = 0.2, nperm = 1, OOB = TRUE, pre1.0_0 = conditional)
{
  
  response <- object@responses
  if (length(response@variables) == 1 && 
      inherits(response@variables[[1]], "Surv"))
    return(varimpsurv(object, mincriterion, conditional, threshold, nperm, OOB, pre1.0_0))
  input <- object@data@get("input")
  xnames <- colnames(input)
  inp <- initVariableFrame(input, trafo = NULL)
  y <- object@responses@variables[[1]]
  if(length(response@variables) != 1)
    stop("cannot compute variable importance measure for multivariate response")
  
  if (conditional || pre1.0_0) {
    if(!all(complete.cases(inp@variables)))
      stop("cannot compute variable importance measure with missing values")
  }
  CLASS <- all(response@is_nominal)
  ORDERED <- all(response@is_ordinal)
  if (CLASS) {
    error <- function(x, oob) mean((levels(y)[sapply(x, which.max)] != 
        y)[oob])
  }
  else {
    if (ORDERED) {
      error <- function(x, oob) mean((sapply(x, which.max) != 
          y)[oob])
    }
    else {
      error <- function(x, oob) mean((unlist(x) - y)[oob]^2)
    }
  }
  
  w <- object@initweights
  if (max(abs(w - 1)) > sqrt(.Machine$double.eps))
    warning(sQuote("varImp"), " with non-unity weights might give misleading results")
  
  ## list for several permutations
  perror <- matrix(0, nrow = nperm*length(object@ensemble), ncol = length(xnames))
  ## this matrix is initialized with values 0 so that a tree that does not 
  ## contain the current variable adds importance 0 to its average importance
  colnames(perror) <- xnames
  for (b in 1:length(object@ensemble)){
    tree <- object@ensemble[[b]]
    
    
    ## if OOB == TRUE use only oob observations, otherwise use all observations in learning sample
    if(OOB){oob <- object@weights[[b]] == 0} else{ oob <- rep(TRUE, length(y))}
    p <- .R_predict(tree, inp, mincriterion, -1L)
    eoob <- error(p, oob)
    
    ## for all variables (j = 1 ... number of variables) 
    for(j in unique(varIDs(tree))){
      for (per in 1:nperm){
        
        if (conditional || pre1.0_0) {
          tmp <- inp
          ccl <- create_cond_list(conditional, threshold, xnames[j], input)
          if (length(ccl) < 1) {
            perm <- sample(which(oob))
          } else {
            perm <- conditional_perm(ccl, xnames, input, tree, oob)
          }
          tmp@variables[[j]][which(oob)] <- tmp@variables[[j]][perm]
          p <- .R_predict(tree, tmp, mincriterion, -1L)
        } else {
          p <- .R_predict(tree, inp, mincriterion, as.integer(j))
        }
        ## run through all rows of perror
        perror[(per+(b-1)*nperm), j] <- (error(p, oob) - eoob)
        
      } ## end of for (per in 1:nperm)
    } ## end of for(j in unique(varIDs(tree)))
  } ## end of for (b in 1:length(object@ensemble))
  
  perror <- as.data.frame(perror)
  #return(MeanDecreaseAccuracy = perror) ## return the whole matrix (= nperm*ntree values per variable)
  return(MeanDecreaseAccuracy = colMeans(perror)) ## return only averages over permutations and trees
}
