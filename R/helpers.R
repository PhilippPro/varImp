# for the current variable of interest, xname,
# create the list of variables to condition on:

create_cond_list <- function(cond, threshold, xname, input) {
  
  stopifnot(is.logical(cond))
  if (!cond) return(NULL)
  if (threshold > 0 & threshold < 1) {
    ctrl <- ctree_control(teststat = "quad", testtype = "Univariate", stump = TRUE)
    xnames <- names(input)
    xnames <- xnames[xnames != xname]
    ct <- ctree(as.formula(paste(xname, "~", paste(xnames, collapse = "+"), collapse = "")),
      data = input, controls = ctrl)
    crit <- ct@tree$criterion[[2]]
    crit[which(is.na(crit))] <- 0
    return(xnames[crit > threshold])
  }
  stop()
}

### extract ID of _all_ variables the tree uses for splitting
varIDs <- function(node) {
  
  v <- c()
  foo <- function(node) {
    if (node[[4]]) return(NULL)
    v <<- c(v, node[[5]][[1]])
    foo(node[[8]])
    foo(node[[9]])
  }
  foo(node)
  return(v)
}
