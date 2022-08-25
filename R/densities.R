## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

LogAdd <- function(X) {
     ##  Calculates log( sum(exp(x)) )  without 'leaving' log space
    if (is.vector(X)) {
        mix <- which.max(X)
        max <- X[mix]
        res <- max + log(sum(exp(X - max)))
    }
    
    if (is.matrix(X)) {
        mv <- apply(X, 1, max)
        res <- mv + log(rowSums(exp(X - mv)))
    }
    
    return(res)
}

dDirichlet <- function(X, alpha, log.p = FALSE) {
    if (isTRUE(all.equal(sum(X), 1))) {
        ix <- is.finite(log(X))
        X[!ix] <- .Machine[["double.xmin"]]
        
        log.lik <- sum(((alpha - 1) * log(X)))
        log.z <- sum(lgamma(alpha)) - lgamma(sum(alpha))
        log.dens <- sum(log.lik) - log.z
    } else {
        log.dens <- -Inf
    }
    
    if (log.p) {
        return(log.dens)
    } else {
        return(exp(log.dens))
    }
}
