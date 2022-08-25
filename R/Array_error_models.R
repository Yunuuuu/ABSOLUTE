## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

AllelicAtten <- function(r, at) {
  return(r * (1 + at) / (1 + at * r))
}

AllelicInvAtten <- function(r, at) {
  return(r / (1 + at) - (at * r))
}

HTx <- function(x, a, b) {
  asinh(a + b * x)
}

HTxInv <- function(x, a, b) {
  return((sinh(x) - a) / b)
}

HTxVar <- function(x, a, b) {
  return(b / (sqrt(1 + (a + b * x)^2)))
}

AllelicTxData <- function(hscn.params, d) {
   tx.b <- (sqrt(exp(hscn.params[["sigma.eta"]]^2) - 1)) / hscn.params[["sigma.nu"]]
   tx.a <- 0
   d.tx <- HTx(d, tx.a, tx.b)

   return(d.tx) 
}


AllelicInvTxData <- function(hscn.params, d) {
   tx.b <- (sqrt(exp(hscn.params[["sigma.eta"]]^2) - 1)) / hscn.params[["sigma.nu"]]
   tx.a <- 0
   d.tx <- HTxInv(d, tx.a, tx.b)

   return(d.tx) 
}


CombAtLL <- function(par, init.par, q, obs, dom1, dom2, lambda.qz.res, pz,
                     sigma.h, verbose=FALSE) {
   if (any(is.na(par))) {
     if (verbose) {
       cat("$")
     }
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
   at <- par[3]

   obs[["error.model"]][["fit.at"]] <- at
   comb <- GetCopyRatioComb(q, delta, b, obs[["error.model"]])

   ll.1 <- CalcNormLoglik(obs[["d.tx"]], obs[["d.stderr"]], obs[["W"]],
                            comb, lambda.qz.res, pz, sigma.h)

   res <- GetAlphaAndTau(b, delta)
   alpha <- res[["alpha"]]
   tau <- res[["tau"]]

   ll.2 <- sum(dnorm(c(alpha, tau, at), init.par, c(0.05, 0.05, 0.01), log=TRUE))

   ll <- ll.1 + ll.2 

   return(-ll)
}

AttenuationOpt <- function(cur.par, obs, dom1, dom2, lambda.qz.res,
                           pz, sigma.h, q, verbose=FALSE) {
   ll <- CombLL(cur.par, q, obs, dom1, dom2, lambda.qz.res,
                 pz, sigma.h )

   if (is.finite(ll)) {
     b <- cur.par[1]
     delta <- exp(cur.par[2])
     res <- GetAlphaAndTau(b, delta)
     alpha <- res[["alpha"]]
     tau <- res[["tau"]]

     ## FIXME: Another AT/at issue
     if ("AT" %in% names(obs$error.model)) {
       obs$error.model$at <- obs$error.model$AT
     }
     at <- obs[["error.model"]][["at"]]
     init.par <- c(alpha, tau, at)

     res  <- optim(par=c(cur.par,at), fn=CombAtLL, gr=NULL,
                   method="Nelder-Mead", q=q, obs=obs, dom1=dom1,
                   dom2=dom2, lambda.qz.res=lambda.qz.res, pz=pz,
                   sigma.h=sigma.h, control=list("maxit"=1000 ),
                   init.par=init.par, verbose=verbose)

      val <- list(b=res[["par"]][1], delta=res[["par"]][2],
                  AT=res[["par"]][3])
     if (verbose) {
       cat("+")
     }
   } else {
     if (verbose) {
       cat("!")
     }
     val <- rep(NA, 3)
   }
   
   return(val)
}

AllelicGetScnaStderrGridDensity <- function(obs, grid, sigma.h, i=NA) {
   sigma.nu <- obs[["error.model"]][["sigma.nu"]]
   sigma.eta <- obs[["error.model"]][["sigma.eta"]]

   b <- (sqrt(exp(sigma.eta^2) - 1)) / sigma.nu
   a <- 0
   N <- length(obs[["d.tx"]])
   grid.dens <- matrix(NA, nrow=N, ncol=length(grid))
   dx <- c(0, diff(grid))

   for (i in seq_len(N)) {
     grid.dens[i, ] <- dnorm(HTx(grid, a, b), obs[["d.tx"]][i],
                             obs[["d.stderr"]][i] + sigma.h) *
                               HTxVar(grid, a, b) * dx      
   }
   
   return(grid.dens)
}
