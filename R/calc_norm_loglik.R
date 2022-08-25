## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

CalcNormLoglik <- function(d, sigma, w, comb, lambda.qz.res, pz, sigma.h, 
                           comb.s=1) {
  if (is.list(lambda.qz.res)) {
    log.pqz <- t((LambdaExpr(lambda.qz.res[["Lambda"]], w))) +
      log(lambda.qz.res[["norm.term"]])
    
    q <- length(lambda.qz.res[["Lambda"]])
    comb <- comb[1:q]
  } else {
    q <- length(comb)
           
    log.pqz <- matrix(log(pz), nrow = length(w), ncol = q + 1, byrow = TRUE)
  }
  
  log.ss.p.1 <- sapply(comb, dnorm, d, sqrt(sigma^2 + comb.s^2 * sigma.h^2),
                       log=TRUE)
  
  ## outlier dist is uniform on [0,3]
  log.comb.seg.mat <- cbind(log.pqz[, c(1:q)] + log.ss.p.1,
                            log(1 / 7) + log.pqz[, q + 1])  
  LL <- sum(LogAdd(log.comb.seg.mat))
  
  return(LL)
}

GetThetaQzPost <- function(d, sigma, w, comb, sigma.h, pi_theta_qz, theta.qz = NA, comb.s = 1) {
    q <- length(comb)
    unif.qz <- rep(1 / (q + 1), q + 1)
    
    log.ss.p.1 <- sapply(comb, dnorm, d,
                         sqrt(sigma^2 + comb.s^2 * sigma.h^2), log = TRUE)
    
    if (is.na(theta.qz)) {
        log.pqz <- matrix(log(unif.qz), nrow = length(w), ncol = q + 1, byrow = TRUE)
    } else {
        log.pqz <- matrix(log(theta.qz), nrow = length(w), ncol = q + 1, byrow = TRUE)
    }
    
    theta.q <- unif.qz[c(1:q)]
    theta.q <- theta.q / sum(theta.q)
##    log_pQ <- matrix(log(theta.q), nrow = length(w), ncol = q, byrow = TRUE)
    
    log.comb.seg.mat <- cbind(log.pqz[, c(1:q)] + log.ss.p.1, log(1/7) +
                              log.pqz[, q + 1])  
    seg_qz <- exp(log.comb.seg.mat - LogAdd(log.comb.seg.mat))    
  
    prior_w = pi_theta_qz$PC
    pi_w = pi_theta_qz$W
    pi_w = pi_theta_qz$W / sum(pi_w)
    
    post_theta_qz = (1 - prior_w) * colSums(w * seg_qz) + (prior_w * pi_w)
    
    return(post_theta_qz)
}

GetSegQzPost <- function(d, sigma, w, comb, lambda.qz.res, sigma.h, 
                         theta.qz=NA, comb.s=1) {
  if (is.list(lambda.qz.res)) {
    log.pqz <- t((LambdaExpr(lambda.qz.res[["Lambda"]], w))) +
      log(lambda.qz.res[["norm.term"]])
    
    ## length of lambda is N_states - 1 ;  == Q
    q <- length(lambda.qz.res[["Lambda"]])
    comb <- comb[1:q]
    log.pq <- log.pqz[, c(1:q)]
    log.pq <- log.pq - LogAdd(log.pq)  ## renormalize
  } else {
    q <- length(comb)   
    theta.q <- theta.qz[c(1:q)]
    theta.q <- theta.q / sum(theta.q)
    log.pqz <- matrix(log(theta.qz), nrow = length(w), ncol = q + 1, byrow = TRUE)
    log.pq <- matrix(log(theta.q), nrow = length(w), ncol = q, byrow = TRUE)
  }
  
  log.ss.p.1 <- sapply(comb, dnorm, d,
                       sqrt(sigma^2 + comb.s^2 * sigma.h^2), log = TRUE)
  
  log.comb.seg.mat <- cbind(log.pqz[, c(1:q)] + log.ss.p.1,
                            log(1/7) + log.pqz[, q + 1])

  seg_qz <- exp(log.comb.seg.mat - LogAdd(log.comb.seg.mat))
  
  log.comb.seg.mat <- log.pq + log.ss.p.1  ## Nseg x Q
  seg.Q <- exp(log.comb.seg.mat - LogAdd(log.comb.seg.mat))
  
  return(list(QZ = seg_qz, Q = seg.Q))
}
