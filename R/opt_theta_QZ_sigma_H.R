OptThetaQzSigmaH <- function(obs, comb, sigma.h.dom, pi_theta_qz, verbose=FALSE) {
  objective <- function(par, obs, comb, lambda) {
    sigma.h <- par
    max.q <- length(comb)
    theta.qz.hat <- GetThetaQzPost(obs[["d.tx"]], obs[["d.stderr"]],
                                   obs[["W"]], comb, sigma.h, 
                                   pi_theta_qz)
    theta.qz.map <- theta.qz.hat
    theta.qz.map <- theta.qz.map / sum(theta.qz.map)
    theta.qz.hat <- theta.qz.map
    if (verbose) {
      cat("S")
    }
    LL <- CalcNormLoglik(obs[["d.tx"]], obs[["d.stderr"]], obs[["W"]], comb, NA, 
                         theta.qz.hat, sigma.h)
    if (is.nan(LL)) {
      LL <- -Inf
    }
    
    return(LL)
  }
  
  max_q = length(comb)
  ## find a value of Lambda to initialize with
  use.sigma.h <- mean(sigma.h.dom)
  theta.qz.hat <- GetThetaQzPost(obs[["d.tx"]], obs[["d.stderr"]],
                                 obs[["W"]], comb, use.sigma.h, pi_theta_qz)
  lambda.qz.res <- GetLambda(theta.qz.hat, obs[["W"]], max.q + 1,
                             verbose=verbose)
  init.lambda <- lambda.qz.res[["Lambda"]]
  
  ## now optimize sigma.h
  res <- optimize(f = objective, interval = sigma.h.dom, tol = 0.001, maximum = TRUE, 
                  obs = obs, comb = comb, lambda = init.lambda)
  
  sigma.h.hat <- res[["maximum"]]
  
  ## now get LL, Lambda, and theta.qz.hat conditional on optimized sigma.h
  theta.qz.hat <- GetThetaQzPost(obs[["d.tx"]], obs[["d.stderr"]],
                                 obs[["W"]], comb, sigma.h.hat, pi_theta_qz)
  
  lambda.qz.res <- GetLambda(theta.qz.hat, obs[["W"]], max.q + 1, verbose=verbose)
  LL <- CalcNormLoglik(obs[["d.tx"]], obs[["d.stderr"]], obs[["W"]],
                       comb, lambda.qz.res, NA, sigma.h.hat)
  
  if (is.nan(LL)) {
    LL <- -Inf
    if (verbose) {
      print("Warning: NaN loglik in opt_theta_Z_sigma_H")
    }
  }

  if (verbose) {
    cat(paste("	sigma.h.hat=", round(sigma.h.hat, 5), sep = ""))
    cat("\n")
  }
  
  return(list(LL = LL, sigma.h.hat = sigma.h.hat, theta_qz_hat = theta.qz.hat, 
              lambda.qz.res = lambda.qz.res))
}
