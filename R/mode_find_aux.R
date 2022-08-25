## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

FindLocationModes <- function(obs, mut, Q, theta.qz, sigma.h, 
                              tau.dom, verbose=FALSE) {
  kAlphaDom <- c(0, 1)
  kDom1 <- c(0, 1)
  kDom2 <- log(c(0.08, 1.05))
  
  ## lambda results
  lambda.qz.res <- GetLambda(theta.qz, obs[["W"]], Q + 1, verbose=verbose)
  pz <- theta.qz[Q + 1]
  
  mode.tab <- MargModeFinder(obs, mut, kDom1, kDom2, Q, lambda.qz.res, pz, sigma.h,
                             tau.dom, verbose=verbose)
  
  if (nrow(mode.tab) == 0) {
    return(list(mode.flag = "DELTA_B_DOM"))
  }
  
  mode.tab <- cbind(mode.tab, NA, NA, NA)
  
  for (i in 1:nrow(mode.tab)) {
    par <- mode.tab[i, c(1, 2)]
    if (obs[["platform"]] == "ARRAY") {
      cur.par = par
      if (is.null(obs$error.model$at)) {
         ## We should not ever be here, if found ,trace it back to source
	 stop("at/AT issue, contact Jeff Gentry jgentry@broadinstitute.org")
      }
      mode.params = list(b=par[1], delta=par[2], AT=obs$error.model$at)
      obs$error.model$fit.at = obs$error.model$at
    } else {
      cur.par <- par
      mode.params <- list(b = par[1], delta = par[2], AT = NA)
    }
    
    LL <- CombLL(cur.par, Q = Q, obs = obs, dom1 = kDom1, dom2 = kDom2,
                  lambda.qz.res, pz, sigma.h)
    
    mode.hess <- hessian(CombLL, x = cur.par, method = "Richardson", Q = Q, 
                         obs = obs, dom1 = kDom1, dom2 = kDom2,
                         lambda.qz.res = lambda.qz.res, pz = pz, 
                         sigma.h = sigma.h)
    
    if (!is.na(mode.hess[1])) {
      mode.curv <- CalcModeLogCurv(par, mode.hess, verbose=verbose)
    } else {
      LL <- NA
      mode.curv <- NA
    }
    
    mode.tab[i, ] <- c(mode.params[["b"]], mode.params[["delta"]],
                       mode.params[["AT"]], LL, mode.curv)
  }
  
  b <- mode.tab[, 1]
  delta <- exp(mode.tab[, 2])
  at <- mode.tab[, 3]
  res <- GetAlphaAndTau(b, delta)
  alpha <- res[["alpha"]]
  tau <- res[["tau"]]
  
  LL <- mode.tab[, 4]
  mode.curv <- mode.tab[, 5]
  mode.tab <- cbind(alpha, tau, at, b, delta, LL, mode.curv)
  colnames(mode.tab) <- c("alpha", "tau", "AT", "b", "delta", "LL", "mode_curv")

  if (verbose) {
    print(mode.tab)
  }
  
  mode.ix <- (mode.tab[, "alpha"] >= kAlphaDom[1] &
              mode.tab[, "alpha"] <= kAlphaDom[2] & 
              mode.tab[, "tau"] >= tau.dom[1] &
              mode.tab[, "tau"] <= tau.dom[2])

  if (verbose) {
    print(paste("removing ", sum(!mode.ix), " / ",
                length(mode.ix), " modes outside of alpha/tau range.", 
                sep = ""))
  }
  mode.tab <- mode.tab[mode.ix, , drop = FALSE]
  
  if (nrow(mode.tab) == 0) {
    return(list(mode.flag = "ALPHA_TAU_DOM"))
  }
  
  return(list(mode.tab = mode.tab))
}

MargModeFinder <- function(obs, mut, dom1, dom2, Q, lambda.qz.res, pz, sigma.h, 
                           tau_dom, b.res=0.125, d.res=0.125, verbose=FALSE) {
  b.grid <- seq(dom1[1], dom1[2], b.res)
  d.grid <- seq(dom2[1], dom2[2], d.res)
  n.b <- length(b.grid)
  n.d <- length(d.grid)
  
  mode.tab <- array(NA, dim = c(n.b * n.d, 3))
  for (i in seq_len(n.b)) {
    for (j in seq_len(n.d)) {
      cur.par <- c(b.grid[i], d.grid[j])
      res <- RunOpt(cur.par, obs, dom1, dom2, lambda.qz.res,
                    pz, sigma.h, Q, verbose=verbose)
      if (!is.na(res)) {
        mode.tab[(i - 1) * n.d + j, ] <- c(res[[1]], res[[2]], res[[3]])
      }
    }
    if (verbose) {
      cat("\n")
    }
  }
  
  ## Try 1d opt for pure tumors
  delta_dom = log(c(1 / tau_dom[2] - 0.05, 1))
  res_1d = run_1d_opt(obs, delta_dom, d.res, lambda.qz.res, pz, sigma.h, Q, verbose=verbose)  
  if (verbose) {
    print("1d mode opt: ")
    print(res_1d)
  }
  mode.tab = rbind(mode.tab, res_1d)
    
  # try opt on SNVs only for diploid tumors
  if (!is.na(mut)) {
    alpha_dom = c(0.1, 1) 
    res_snv_only = run_diploid_snv_purity_opt(obs, mut$mut_CN_dat, 
                                              mut$Pi_som_theta_Q, mut$mut_class_W, 
                                              alpha_dom, verbose=verbose)
    
    if (!is.na(res_snv_only)) {
      if (verbose) {
        print("SNV_only opt: ")
        print(res_snv_only)
      }
      mode.tab = rbind(mode.tab, res_snv_only)
    }
  }
    
  ## find unique modes in table
  ix <- !is.na(mode.tab[, 1])
  mode.list <- mode.tab[ix, c(1, 2), drop = FALSE]
  umodes = unique(round(mode.list, 2))
  
  if (verbose) {
    print(paste(nrow(umodes), " unique modes found", sep = ""))
  }
  umodes <- umodes[umodes[, 1] >= dom1[1] & umodes[, 2] >= dom2[1], , drop = FALSE]
  if (verbose) {
    print(paste(nrow(umodes), " modes in b / delta range.", sep = ""))
  }
  return(umodes)
}

GetUniqModes <- function(mode.list, tol) {
    N <- nrow(mode.list)
    umat <- matrix(FALSE, nrow = N, ncol = N)
    
    if (N < 2) {
      return(mode.list)
    }
    
    for (i in 2:N) {
      for (j in 1:(i - 1)) {
        umat[i, j] <- isTRUE(all.equal(mode.list[i, ], mode.list[j, ], tol = tol))
      }
    }
    
    uniq.ix <- rowSums(umat) == 0
    
    return(mode.list[uniq.ix, , drop = FALSE])
}


RunOpt <- function(cur.par, obs, dom1, dom2, lambda.qz.res, pz, sigma.h, 
                   Q, eval.hessian=FALSE, verbose=FALSE) {
  LL <- CombLL(cur.par, Q = Q, obs = obs, dom1 = dom1, dom2 = dom2,
               lambda.qz.res, pz, sigma.h)
    
  if (is.finite(LL)) {
    res <- optim(par = cur.par, fn = CombLL, gr = NULL, method = "Nelder-Mead", 
                 Q = Q, obs = obs, dom1 = dom1, dom2 = dom2,
                 lambda.qz.res = lambda.qz.res, 
                 pz = pz, sigma.h = sigma.h, control = list(maxit = 1000),
                 hessian = eval.hessian)
    res <- list(res[["par"]][1], res[["par"]][2], -res[["value"]],
                res[["hessian"]])
    if (verbose) {
      cat(".")
    }
  } else {
    if (verbose) {
      cat("!")
    }
    res <- NA
  }
  
  return(res)
}

run_1d_opt = function(obs, dom, d_res, lambda_qz_res, pz, sigma_h, Q, verbose=FALSE) {
  comb_1d_ll = function(par, Q, obs, lambda_qz_res, pz, sigma_h) {
    delta = exp(par)
    b = 0
    
    comb =  GetCopyRatioComb(Q, delta, b, obs$error.model)
    LL = CalcNormLoglik(obs$d.tx, obs$d.stderr, obs$W, comb, lambda_qz_res, pz, sigma_h)
    
    return(-LL)
  }
  
  d_grid = seq(dom[1], dom[2], d_res)
  mode_tab = array(NA, dim=c( length(d_grid), 3))
  
  for (i in seq_along(d_grid)) {
    opt = nlm(f=comb_1d_ll, p=d_grid[i], Q=Q, obs=obs, lambda_qz_res=lambda_qz_res, pz=pz, sigma_h=sigma_h)
    
    mode_tab[i, ] = c(0, opt$estimate, -opt$minimum)
    if (verbose) {
      cat("-")
    }
  }

  if (verbose) {
    cat("\n")
  }
  
  return(mode_tab)
}


CombLL <- function(par, Q, obs, dom1, dom2, lambda.qz.res, pz, sigma.h) {
    if (any(is.na(par))) {
        cat("$")
        return(Inf)
    }
    
    if (par[1] < -0.05 | par[1] > dom1[2]) {
        return(Inf)
    }
    if (par[2] < log(0.07) | par[2] > dom2[2]) {
        return(Inf)
    }
    
    b <- par[1]
    delta <- exp(par[2])
    
    comb <- GetCopyRatioComb(Q, delta, b, obs[["error.model"]])
    LL <- CalcNormLoglik(obs[["d.tx"]], obs[["d.stderr"]], obs[["W"]], comb,
                         lambda.qz.res, pz, sigma.h)
    
    return(-LL)
}
