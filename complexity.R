require(zoo) # for rollapply
require(dplyr) #
require(scales) # rescale
require(psych) # MSSD

complexity <- function(x, scaleMin = min(x, na.rm = T), scaleMax = max(x, na.rm = T), width = 7, measure = "complexity", rescale = FALSE) {

if(!is.numeric(x)){return("Please provide a numeric vector.")}

# MSSD function taken from the psych library ----


# RMSSD of the segment, standardized by the maximum possible RMSSD ----
  fluctDegree <- function(x, scaleMin, scaleMax) {

    seq <- length(x)
    fmax <- mssd(rep_len(c(scaleMin, scaleMax), length.out = seq)) # Maximum possible MSSD. Vor (very) long time series, this would be scaleMax - scaleMin, but for smaller windows, it is larger.
    fobs <- mssd(x, lag = 1)
    f <- fobs / fmax
    return(f)
    }

# Deviation from a hypothetical uniform distribution
  distDegree <- function(x, scaleMin, scaleMax) {
    uniform <- seq(from = scaleMin, to = scaleMax, length.out = length(x)) # Make a hypothetical distribution for the window
    empirical <- sort(x) # Rank-order observed values
    uni.diff <- diff(uniform) # Calculate differences in hypothetical distribution
    emp.diff <- diff(empirical) # Calculate differences in observed distribution
    deviation <- uni.diff - emp.diff# Calculate deviation
    dev.h <- deviation * ((sign(deviation) + 1)/2) # Apply a Heaviside step function to eliminate negatives
    div.diff <- dev.h / uni.diff
    #print(normed) # debug
    D <- 1 - mean(div.diff) # If there were no deviations from the uniform distribution, this would be 0
    return(D)}

  fluctuation <- rollapply(x, width = width, FUN = fluctDegree, scaleMin = scaleMin, scaleMax = scaleMax, partial = F, fill = NA)
  distribution <- rollapply(x, width = width, FUN = distDegree, scaleMin = scaleMin, scaleMax = scaleMax, partial = F, fill = NA)
  complexity <- fluctuation * distribution

  if(rescale){
    if(measure == "distribution"){return(scales::rescale(distribution, to = c(scaleMin, scaleMax)))}
    if(measure == "fluctuation"){return(scales::rescale(fluctuation, to = c(scaleMin, scaleMax)))}
    if(measure == "complexity"){return(scales::rescale(complexity, to = c(scaleMin, scaleMax)))}
  }


  if(measure == "distribution"){return(distribution)}
  if(measure == "fluctuation"){return(fluctuation)}
  if(measure == "complexity"){return(complexity)}
}
