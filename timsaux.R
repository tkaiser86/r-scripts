# Tim's miscellaneous R functions
# Included so far:
# - Function for converting among effect sizes (r, d, LogOdds and from these to BESD)
# - omega-squared for aov objects from https://stats.stackexchange.com/questions/2962/omega-squared-for-measure-of-effect-in-r

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
    res = .50 + abs(x)/2
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
