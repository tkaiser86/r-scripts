# Simple miscellaneous R functions. Some of those might be included in other packages. In this case, I was not aware of that. 
# Included so far:
# - omega-squared for aov objects from https://stats.stackexchange.com/questions/2962/omega-squared-for-measure-of-effect-in-r
# - Mean absolute difference (mad), an alternative to standard deviation
# - Distance correlation for a matrix
# - Probability of obtaining two statistically different test scores, classical test theory version (PDTSctt)
#
# source("https://raw.githubusercontent.com/tkasci/r-scripts/master/timsaux.R")

convert.effect <- function(x, from, to, a = 4){
  if(from == "logodds" & to == "d"){
    res <- x * sqrt(3)/pi
    return(res)
  }
  if(from == "d" & to == "logodds"){
    res <- x * pi/sqrt(3)
    return(res)
  }
  if(from == "r" & to == "d"){
    res <- 2*x / sqrt(1-x^2)
    return(res)
  }
  if(from == "d" & to == "r"){
    res <- x / sqrt(x^2+a)
    return(res)
  }
  if(from == "r" & to == "besd"){
    res = .50 + x/2
    return(res)
  }
  if(from == "d" & to == "besd"){
    y <- x / sqrt(x^2+a)
    res <- .50 + y/2
    return(res)
  }
  if(from == "logodds" & to == "besd"){
    y <- x * sqrt(3)/pi  # logodds to d
    y <- y / sqrt(y^2+a) # d to r
    res <- .50 - y/2
    return(res)
  }
}

omega_sq <- function(aov_in, neg2zero=T){
  aovtab <- summary(aov_in)[[1]]
  n_terms <- length(aovtab[["Sum Sq"]]) - 1
  output <- rep(-1, n_terms)
  SSr <- aovtab[["Sum Sq"]][n_terms + 1]
  MSr <- aovtab[["Mean Sq"]][n_terms + 1]
  SSt <- sum(aovtab[["Sum Sq"]])
  for(i in 1:n_terms){
    SSm <- aovtab[["Sum Sq"]][i]
    DFm <- aovtab[["Df"]][i]
    output[i] <- (SSm-DFm*MSr)/(SSt+MSr)
    if(neg2zero & output[i] < 0){output[i] <- 0}
  }
  names(output) <- rownames(aovtab)[1:n_terms]
  
  return(output)
}

# Average absolute deviation, an alternative to the standard deviation
aad <- function(x){
  sum(abs(x-mean(x))) / length(x)
}

# Taleb-style correction of a correlation coefficient
rescale.cor <- function(rho){
  abs(rho) * (1-(sqrt(1-rho^2)))
}

# distance correlation (?energy::dcor) applied to a data matrix or frame
dcor.mat <-
  function(data,
           bias.correct = F,
           square.root = T,
           ignore.warning = F) {
    require(energy)
    if (ignore.warning == F) {
      if (nrow(data) > 3000) {
        return(
          "Distance correlation matrices for data sets longer than 3000 is memory-intensive and calculation time rises exponentially. If you really want to do this, set ignore.warning = T."
        )
      }
    }
    results <- matrix(nrow = ncol(data), ncol = ncol(data))
    colnames(results) <- colnames(data)
    rownames(results) <- colnames(data)
    
    for (i in 1:ncol(results)) {
      for (j in 1:nrow(results)) {
        if (bias.correct) {
          results[i, j] <- bcdcor(x = data[, i], y = data[, j])
        } else {
          results[i, j] <- dcor(x = data[, i], y = data[, j])
        }
        
      }
    }
    if (square.root & bias.correct) {
      results[which(results < 0)] <- 0
      results %>% sqrt() -> results
    }
    results
  }

# fast distance correlation (https://arxiv.org/pdf/1810.11332v1.pdf)
fastDcov <- function(x, y){
  n <- length(x)
  x <- sort(x)
  
  si <- cumsum(x)
  s <- si[n]
  a.x <- (-(n-2)/2/n) * x + (s-2*si)
  a.x
}

pdts_ctt <- function(x, rel = NULL, alpha=.975){
  if(is.null(rel)){stop("Please specify the reliability coefficient.")}
  zscore <- qnorm(alpha)
  LSD <- zscore * sd(x, na.rm=T) * sqrt(2*(1-rel))
  score.combinations <- combn(x, m=2)
  diff.scores <- abs(score.combinations[1,] - score.combinations[2,])
  diff.tab <- table(diff.scores > LSD) 
  diff.tab[2] / sum(diff.tab)
}

maximize.pdts <-  function(x, length, trials = 1000, verbose=F, reliability = NULL){
  pdts.selections <- as.data.frame(matrix(nrow=trials, ncol=(length+3)))
  
  if(length > ncol(x)){stop("Length of subscale is larger than the scale itself!")}
  indices <- 1:ncol(x)
  for(i in 1:trials){
    if(verbose){message(i)}
    subselection <- sample(indices, length, replace = F)
    pdts.selections[i,1:length] <- subselection
    subset.data <- x[,subselection]
    
    alpha.val <- psych::alpha(subset.data)$total$raw_alpha
    omega.val <- psych::omega(subset.data, plot=F)$omega_h
    if(is.null(reliability)){rel <- omega.val} else {rel <- reliability}
    pdts.selections[i,(length+1)] <- pdts_ctt(rowMeans(subset.data), rel = rel)
    pdts.selections[i,(length+2)] <- alpha.val
    pdts.selections[i,(length+3)] <- omega.val}
  
  
  return(pdts.selections)
}
