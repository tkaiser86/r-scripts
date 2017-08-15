# Convert a lag 1 VAR model from vars to edge lists processable by qgraph. Bonferroni-corrects all coefficients.
# nonsig: boolean. Include non-significant parameters in edge lists.
# alpha: alpha level 
# bonferroni: Bonferroni-corrects alpha levels for temporal and contemporal associations
# granger: include only Granger-causal temporal associations
require(vars)
VARtoEdges <- function(x, nonsig = F, alpha = .05, bonferroni = T, granger = F){
  VARsum <- summary(x)
  
  if(bonferroni){
    alphacor <- alpha / ((x$K*(x$K-1))/2) # Correct for every pair
    alphareg <- alpha / x$K # Correct for every predictor
  } else {
    alphacor = alpha
    alphareg = alpha
  }
  
  edgelist <- list()
  temporal.edges <- matrix(nrow = x$K, ncol = x$K)
  
  if(!nonsig){significance <- matrix(nrow = x$K, ncol = x$K)}
  for(i in 1:x$K){
    for(j in 1:x$K){
      temporal.edges[j,i] <- VARsum$varresult[[i]]$coefficients[j]
      if(!nonsig){
        significance[j,i] <- VARsum$varresult[[i]]$coefficients[j,4]
      }
    }}
  
  if(!nonsig){
    significance <- significance < alphareg
    temporal.edges <- temporal.edges * significance}
  df.edges <- as.data.frame(temporal.edges)
  rownames(df.edges) <- colnames(x$datamat)[1:x$K]
  colnames(df.edges) <- colnames(x$datamat)[1:x$K]
  edgelist$temporal <- df.edges
  
  contemp.edges <- VARsum$corres
  if(!nonsig){
    significance <- matrix(nrow= x$K, ncol = x$K)
    for(i in 1:x$K){
      for(j in 1:x$K){
        correlation <- VARsum$corres[j,i]
        tVal <- correlation / sqrt((1-correlation^2) / (x$obs-2))
        pVal <- 1 - pt(tVal, x$obs-2)
        significance[j,i] <- pVal < alphacor
      }}
    contemp.edges <- contemp.edges * significance
  }
  
  if(granger){
    granger <- NULL
    varnames <- colnames(VARsum$corres)
    for(i in 1:x$K){
      granger[i] <- vars::causality(x, cause = varnames[i])$Granger$p.value
    }
    granger.causality <- matrix(rep(granger, ols.var$K), ncol = 5) < .05
    contemp.edges <- contemp.edges * granger.causality
  }
  edgelist$contemp <- contemp.edges
  return(edgelist)
}
