## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

ClassifySamplesWgdByProfile <- function(segobj) {
    CalcHighM <- function(w, seg.ix, seg.qz) {
        uniq.seg.ix <- unique(seg.ix)
        n.seg <- length(uniq.seg.ix)
        high.tab <- array(NA, dim = c(n.seg, ncol(seg.qz)))
        low.tab <- array(NA, dim = c(n.seg, ncol(seg.qz)))       
        found.pairs <- rep(FALSE, n.seg)
        use.w <- rep(NA, n.seg)
        for (i in 1:n.seg) {
            if (sum(seg.ix == uniq.seg.ix[i]) != 2) {
                next
            }
            
            found.pairs[i] <- TRUE
            pair.ix <- which(seg.ix == uniq.seg.ix[i])
            
            use.w[i] <- w[pair.ix[1]]
            high.tab[i, ] <- seg.qz[pair.ix[2], ]
            low.tab[i, ] <- seg.qz[pair.ix[1], ]
        }
        
        high.tab <- high.tab[found.pairs, ]
        use.w <- use.w[found.pairs]
        use.w <- use.w / sum(use.w)
        
        low.m <- colSums(use.w * low.tab)
        high.m <- colSums(use.w * high.tab)
        
        return(rbind(low.m, high.m))
    }
    
    q <- 15
    wgd <- 0
    w <- segobj[["as.seg.dat"]][, "W"]
    seg.ix <- segobj[["as.seg.dat"]][, "seg.ix"]
    seg.q <- segobj[["mode.res"]][["seg.q.tab"]][1, , ]
    
    g.m <- CalcHighM(w, seg.ix, seg.q)
    low.m <- g.m[1, ]
    high.m <- g.m[2, ]
    
    low_p <- sum(low.m * (c(1:q) - 1))
    high_p <- sum(high.m * (c(1:q) - 1))
    
    if (high.m[3] > 0.5) {
        wgd <- 1
    }
    if (sum(high.m[c(3:length(high.m))]) > 0.5) {
        wgd <- 1
    }
    if (sum(high.m[c(4:length(high.m))]) > 0.4) {
        wgd <- 2
    }
    
    return(wgd)
}
