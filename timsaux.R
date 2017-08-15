# Tim's miscellaneous R functions
# Included so far:
# - Function for converting among effect sizes (r, d, LogOdds and from these to BESD)

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

