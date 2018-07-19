
# create_cond_list
#
# create the list of variables to condition on for the current variable of interest, xname,
create_cond_list = function(cond, threshold, xname, input) {
  stopifnot(is.logical(cond))
  if (!cond) return(NULL)
  if (threshold > 0 & threshold < 1) {
    ctrl = ctree_control(teststat = "quad", testtype = "Univariate", stump = TRUE)
    xnames = names(input)
    xnames = xnames[xnames != xname]
    ct = ctree(as.formula(paste(xname, "~", paste(xnames, collapse = "+"), collapse = "")),
      data = input, controls = ctrl)
    crit = ct@tree$criterion[[2]]
    crit[which(is.na(crit))] = 0
    return(xnames[crit > threshold])
  }
  stop()
}

### extract ID of _all_ variables the tree uses for splitting
varIDs = function(node) {
  
  v = c()
  foo = function(node) {
    if (node[[4]]) return(NULL)
    v <<- c(v, node[[5]][[1]])
    foo(node[[8]])
    foo(node[[9]])
  }
  foo(node)
  return(v)
}

conditional_perm = function(cond, xnames, input, tree, oob) {
  
  ## intitial partitioning => all observations in one partition
  parts = rep(1, length(oob))
  
  ## develop partitioning by going over all the conditiong variables
  for (condVar in cond) {
    
    ## varID is variable index or column number of input (predictor matrix) 
    ## not variable name!
    varID = which(xnames == condVar)
    
    ## if conditioning variable is not used for splitting in current tree
    ## proceed with next conditioning variable
    cl = cutpoints_list(tree, varID)
    if (is.null(cl)) next
    
    ## proceed cutpoints for different types of variables
    x = input[, varID]
    xclass = class(x)[1]
    if (xclass == "integer") xclass = "numeric"
    
    block = switch(xclass, "numeric" = cut(x, breaks = c(-Inf, sort(unique(cl)), Inf)),
      "ordered" = cut(as.numeric(x), breaks =  c(-Inf, sort(unique(cl)), Inf)),
      "factor" = {
        CL = matrix(as.logical(cl), nrow = nlevels(x))                            
        rs = rowSums(CL)
        dlev = (1:nrow(CL))[rs %in% rs[duplicated(rs)]]
        fuse = c()
        for (ii in dlev) {
          for (j in dlev[dlev > ii]) {
            if (all(CL[ii,] == CL[j,])) fuse = rbind(fuse, c(ii, j))
          }
        }
        xlev = 1:nlevels(x)
        newl = nlevels(x) + 1
        block = as.integer(x)
        for (l in xlev) {
          if (NROW(fuse) == 0) break
          if (any(fuse[, 1] == l)) {
            f = c(l, fuse[fuse[, 1] == l, 2])
            fuse = fuse[!fuse[,1] %in% f, , drop = FALSE]
            block[block %in% f] = newl
            newl = newl + 1
          }
        }
        as.factor(block)
      })
    ## add partitioning based on the split points the variable to the 
    ## current partitioning
    parts = interaction(parts, as.numeric(block), drop = TRUE, sep = "")
  }
  
  ## if none of the conditioning variables are used in the tree
  if (!length(levels(parts)) > 1) {
    perm = sample(which(oob))
    return(perm)
  } else {
    ## one conditional permutation
    perm = 1:nrow(input)
    for(part in levels(parts)){
      index = which(parts == part & oob)
      if (length(index) > 1)
        perm[index] = index[sample.int(length(index))]
    }
    return(perm[oob])
  }
}

# cutpoints_list() returns:
# - vector of cutpoints (length=number of cutpoints) 
#   if variable is continuous
# - vector of indicators (length=number of categories x number of cutpoints)
#   if variable is categorical (nominal or ordered)
cutpoints_list = function(tree, variableID) {
  
  cutp = function(node) {
    if (node[[4]]) return(NULL)
    cp = NULL
    if (node[[5]][[1]] == variableID)
      cp = node[[5]][[3]]
    nl = cutp(node[[8]])
    nr = cutp(node[[9]])
    return(c(cp, nl, nr))
  }
  return(cutp(tree))
}
