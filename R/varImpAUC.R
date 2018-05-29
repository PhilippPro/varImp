# Berechnung des permutation VIM basierend auf dem AUC für mehrkategoriale Zielgrößen

#' varImpAUC
#'
#' @param object 
#' @param mincriterion 
#' @param conditional 
#' @param threshold 
#' @param nperm 
#' @param OOB 
#' @param pre1.0_0 
#' @param method 
#'
#' @return vector with computed permutation importance for each variable
#' @export
#'
#' @examples 
#' 
#' # multiclass case
#' data(iris)
#' iris.cf <- cforest(Species ~ ., data = iris, control = cforest_unbiased(mtry = 2, ntree = 50))
#' set.seed(123)
#' a = varImpAUC(object = iris.cf, method = "ovo")
#' set.seed(123)
#' b = varImpAUC(object = iris.cf, method = "ova") 
varImpAUC <- function (object, mincriterion = 0, conditional = FALSE, threshold = 0.2, 
  nperm = 1, OOB = TRUE, pre1.0_0 = conditional, method=c("ovo","ova")) { 
  # vgl. Janitza
  method=match.arg(method)
  response <- object@responses
  input <- object@data@get("input")
  xnames <- colnames(input)
  inp <- initVariableFrame(input, trafo = NULL)
  y <- object@responses@variables[[1]]
  if (length(response@variables) != 1) 
    stop("cannot compute variable importance measure for multivariate response")
  if (conditional || pre1.0_0) {
    if (!all(complete.cases(inp@variables))) 
      stop("cannot compute variable importance measure with missing values")
  }
  CLASS <- all(response@is_nominal)
  ORDERED <- all(response@is_ordinal)
  if (CLASS) {
    if (nlevels(y) > 2) {
      if(method=="ova"){ ########################################################### one-versus-all Verfahren 
        error <- function(x, oob) {
          xoob <- t(sapply(x, function(x) x))[oob,]
          yoob <- y[oob]
          return(measures::multiclass.AUNU(xoob, yoob))
        }
      } else if(method=="ovo"){ ############################# one-versus-one, paarweises Verfahren (Hand & Till)
        error <- function(x, oob) {
          xoob <- t(sapply(x, function(x) x))[oob,]
          yoob <- y[oob]
          return(measures::multiclass.AU1U(xoob, yoob))
        }
      }
      ############# AUC-Berechnung für den Fall einer binären Zielgröße (s. Janitza) ############################
    } else { 
      error <- function(x, oob) {
        xoob <- sapply(x, function(x) x[1])[oob]
        yoob <- y[oob]
        pos = levels(y)[1]
        return(measures::AUC(xoob, yoob, positive = pos))
      }
    }
  } else {
    if (ORDERED) {
      error <- function(x, oob) mean((sapply(x, which.max) != y)[oob])
    }
    else {
      error <- function(x, oob) mean((unlist(x) - y)[oob]^2)
    }
  }
  w <- object@initweights
  if (max(abs(w - 1)) > sqrt(.Machine$double.eps)) 
    warning(sQuote("varimp"), " with non-unity weights might give misleading results")
  perror <- matrix(0, nrow = nperm * length(object@ensemble), ncol = length(xnames))
  colnames(perror) <- xnames
  for (b in 1:length(object@ensemble)) {
    tree <- object@ensemble[[b]]
    if (OOB) {
      oob <- object@weights[[b]] == 0
    } else {
      oob <- rep(TRUE, length(xnames))
    }
    p <- .Call("R_predict", tree, inp, mincriterion, -1L, 
      PACKAGE = "party")
    eoob <- error(p, oob)
    for (j in unique(varIDs(tree))) {
      for (per in 1:nperm) {
        if (conditional || pre1.0_0) {
          tmp <- inp
          ccl <- create_cond_list(conditional, threshold, 
            xnames[j], input)
          if (is.null(ccl)) {
            perm <- sample(which(oob))
          }
          else {
            perm <- conditional_perm(ccl, xnames, input, 
              tree, oob)
          }
          tmp@variables[[j]][which(oob)] <- tmp@variables[[j]][perm]
          p <- .Call("R_predict", tree, tmp, mincriterion,-1L, PACKAGE = "party")
        } else {
          p <- .Call("R_predict", tree, inp, mincriterion, 
            as.integer(j), PACKAGE = "party")
        }
        perror[(per + (b - 1) * nperm), j] <- - (error(p,oob) - eoob)
      }
    }
  }
  perror <- as.data.frame(perror)
  return(MeanDecrease = colMeans(perror, na.rm = TRUE))
}
