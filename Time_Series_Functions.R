## Function to group values that are less than or equal to a certain distance from each other

group_ts <- function(x, interval = 1) {
  
  xs <- sort(x)
  gs <- rep(1, length(x))
  gi <- 1
  
  if(length(x) >= 1) {
    
    for (i in 2:length(x)) {
      
      x_diff <- xs[i] - xs[i - 1]
      
      if (x_diff > interval) gi <- gi + 1
      
      gs[i] <- gi
      
    }
    
  }
  
  return(gs[match(x, xs)])
  
}

#' Gaussian temporal kernel smooth with interpolation
#'
#' @param x Input vector to smooth
#' @param t Integer vector of length length(x) corresponding to time steps of observations
#' @param n Number of interpolated time steps in output. If missing, n = length(t)
#' @param bw bandwidth (95% Gaussian quantile range); defaults to 10 time steps
#'
#' @returns A data.frame
gaussian_smooth <- function(x, t, n, bw = 10) {
  
  if (missing(n)) {
    tseq <- t
  } else {
    tseq <- seq(min(t), max(t), length.out = n)
  }
  
  tdiff <- unclass(outer(tseq, tseq, FUN = "-"))
  twght <- exp( -( tdiff^2 / (2 * ((bw/2)/1.96)^2) ) )
  
  xseq <- approx(t, x, xout = tseq)$y
  
  xsmooth <- rep(NA_real_, length(x))
  
  for (i in seq_along(tseq)) {
    
    sub <- which(tseq >= tseq[i] - 20 & tseq <= tseq[i] + 20)
    
    xseq_i <- xseq[sub]
    wght_i <- twght[sub, i] / sum(twght[sub, i])
    
    xsmooth[i] <- weighted.mean(xseq_i, wght_i)
    
  }
  
  return(data.frame(t = tseq, x = xsmooth))
  
}
