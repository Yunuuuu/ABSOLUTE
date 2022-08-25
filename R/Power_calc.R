## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

GetCovForPow <- function(required.power, eps, fdr, delta) {
  max <- 1000
  cov <- 2
  
  if (delta <= 0 | delta > 1) {
    return(NA)
  }
  
  while (TRUE) {
    pow <- PowerCalc(cov, eps, fdr, delta)
    if (pow >= required.power) {
      break
    }
    if (cov > max) {
      cov <- NA
      break
    }
    
    cov <- cov + 1
  }
  return(cov)
}

PowerCalc <- function(n, eps, fdr, delta) {
  p <- dbinom(c(0:n), size = n, prob = eps)
  pv <- 1 - cumsum(c(0, p))
  ks <- min(which(pv <= fdr))
  
  cval <- (fdr - pv[ks]) / (pv[ks - 1] - pv[ks])
  
  p1 <- dbinom(c(0:n), size = n, prob = delta)
  pow <- 1 - sum(p1[c(1:(ks - 1))]) + cval * p1[ks - 1]
  
  return(pow)
}
