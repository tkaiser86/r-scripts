# Convert a lag 1 VAR model from vars to edge lists processable by qgraph. Bonferroni-corrects all coefficients.
# p.adjust = Adjust p-values for VAR models using the methods supported by stats::p.adust
# see ?p.adjust for more information
# granger: include only Granger-causal temporal associations
require(vars)
VARtoEdges <- function(x, p.adjust = "none", granger = F) {
  var.sum <- summary(x)
  nvar <- x$K
  maxp <- x$p
  mat.list <- list()
  
  for (a in 1:maxp) {
    mat.list[[a]] <- matrix(ncol = nvar, nrow = nvar)
    colnames(mat.list[[a]]) <- var.sum$names
    rownames(mat.list[[a]]) <- var.sum$names
    for (b in 1:nvar) {
      coefs <- data.frame(var.sum$varresult[[b]]$coefficients)
      current <- coefs[grep(a, row.names(coefs)), c(1, 4)]
      significant <- p.adjust(current[, 2], method = p.adjust) < .05
      mat.list[[a]][, b] <- current[, 1] * significant
    }
  }
  
  if(granger){
    granger.causality <- vector()
    for(i in 1:length(summary(varmodel)$names)){
      granger.causality[i] <- causality(varmodel, cause = summary(varmodel)$names[i])$Granger$p.value  
    }
    names(granger.causality) <- summary(varmodel)$names
  }
  
  return(mat.list)
}
