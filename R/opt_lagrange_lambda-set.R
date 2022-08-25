## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

L2Loss <- function(lambda, n.states, theta, W) {
  res <- LambdaEval(lambda, W, n.states) 
  loss <- sqrt(sum((theta - res)^2))
  
  if (!is.finite(loss)) {
    stop()
  }  

  return(loss)
}

GetInitGues <- function(W, n.states, theta) {
  l0 <- seq(1, 0, length=(n.states-1))
  N <- length(W)

  l0 <- l0 * runif(1, N, 3*N)
  sgn <- ifelse(rbinom(1, 1, 0.5) == 1, -1, 1)
  return(l0 * sgn)
}

LambdaEval <- function(lambda, W, n.states) {
  log.p.q.terms <- LambdaExpr(lambda, W)  
  p.q.terms <- exp(log.p.q.terms)
  p.q.terms[p.q.terms == Inf] <- .Machine[["double.xmax"]] / 2
  seg.totals <- apply(p.q.terms, 2, sum)
  seg.densities <-  t(p.q.terms) / seg.totals
  g1 <- apply(W * seg.densities, 2, sum) 
  
  if (any(!is.finite(g1))) {
    stop()
  }
  
  return(g1)
}

LambdaExpr <- function(lambda, W) {
  n.states <- length(lambda) + 1
  return(-c(0, lambda) * ((c(1:n.states)-1) %*% matrix(W, nrow=1)))
}

## theta is a vector of length n.states
## need as many lambda's as theta's
CalcLambda <- function(theta, W, n.states, first.lambda.guess,
                       verbose=FALSE) {
  if (all(!is.na(first.lambda.guess))) {
    cur.lambda <- first.lambda.guess
  } else {
    ## init guess
    cur.lambda <-  runif(c(n.states-1), -1e3, 1e3)    
  }
  
  while (!is.finite(L2Loss(cur.lambda, n.states, theta, W))) {
    cur.lambda <- GetInitGues(W, n.states)
    if (verbose) {
      cat("!")
    }
  }

  stol <- 1e-6
  gtol <- 1e-6
  
  res  <- try(nlm(f=L2Loss, p=cur.lambda, n.states=n.states, theta=theta,
                  W=W, steptol=stol, gradtol=gtol), silent=TRUE)
  if (inherits(res, "try-error")) {
    if (verbose) {
      cat("!")
    }
    return(list(resid=Inf))
  }
  
  lambda <- res[["estimate"]]
  resid=res[["minimum"]]
    
  ## calc normalization factor (for each seg)
  log.p.q.terms <- LambdaExpr(lambda, W)
  
  p.q.terms <- exp(log.p.q.terms)
  seg.totals <- apply(p.q.terms, 2, sum)
  lambda.norm.term <- 1 / seg.totals
  
  return(list(Lambda=lambda, norm.term=lambda.norm.term, resid=resid))
}

GetLambda <- function(theta, W, n.states, lambda.init=NA, w.thresh=0,
                      verbose=FALSE) {
  if (all(theta == theta[1])) {
    lambda.res <- list()
    lambda.res[["Lambda"]] <- rep(0, (n.states-1))
    lambda.res[["norm.term"]] <- rep(1 / n.states, length(W))
    lambda.res[["resid"]] <- 0 
    return(lambda.res)
  }

  if (!is.na(lambda.init)) {
    lambda.res <- CalcLambda(theta, W, n.states, first.lambda.guess=lambda.init,
                             verbose=verbose)
    resid <- lambda.res[["resid"]]

    if (resid < 0.05) {
      return(lambda.res)
    }
  }
  wsix <- sort(W, index.return=TRUE)[["ix"]]
  orig.w <- W
  if (w.thresh > 0) {
    n.small <- max(which(cumsum(W[wsix]) < w.thresh))
    small.ix <- wsix[c(1:n.small)]
    
    W <- W[-small.ix]
    W <- W / sum(W)
  }
  
  theta[theta < .Machine$double.xmin] = .Machine$double.xmin
  
  resid <- Inf
  wix <- sort(W, decreasing=TRUE, index.return=TRUE)[["ix"]]
  iix <- 1
  init.ix <- wix[iix]
  
  while (resid > 0.05 & iix <= length(W)) {
    lambda.init <- -1/W[init.ix] * log(theta[c(2:length(theta))] /theta[1]) / c(1:(length(theta)-1))
    lambda.init[lambda.init == Inf] = .Machine$double.xmax / 2
    lambda.init[lambda.init == -Inf] = -.Machine$double.xmax / 2
    
    lambda.res <- CalcLambda(theta, W, n.states, first.lambda.guess=lambda.init,
                             verbose=verbose)
    resid <- lambda.res[["resid"]]
    iix <- iix+1
    init.ix <- wix[iix]
  }

  if (resid > 0.05) {
    if (verbose) {
      print(paste("Warning: Lambda resid = ", resid, sep=""))
    }
    stop()
  }
  
  if (w.thresh > 0) {
    norm.term <- rep(NA, length(orig.w))
    norm.term[small.ix] <- 1/n.states
    norm.term[-small.ix] <- lambda.res[["norm.term"]]
    lambda.res[["norm.term"]] <- norm.term
  }
  
  return(lambda.res)
}
