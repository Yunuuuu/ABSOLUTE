## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

AllelicCalcChrArmDistr <- function(seg.obj, seg.q, chr.arms.dat) {
  n.arm <- nrow(chr.arms.dat)
  chr.arm.tab <- array(NA, dim=c(2, n.arm, ncol(seg.q)))
  
  for (i in 1:n.arm) {
    chr.dat <- GetAllelicChrArmSegs(seg.obj, chr.arms.dat[i, ])      
    
    if (length(chr.dat[["int.w"]]) == 0 ) {
      next
    }
    
    chr.arm <- array(NA, dim=c(2, ncol(seg.q)))
    
    chr.arm[1, ] <- colSums(seg.q[chr.dat[["low.ix"]], , drop=FALSE] *
      chr.dat[["int.w"]])
    chr.arm[2, ] <- colSums(seg.q[chr.dat[["high.ix"]], , drop=FALSE] *
      chr.dat[["int.w"]])
    
    chr.arm.tab[, i, ] <- 0
    chr.arm.tab[1, i, which.max(chr.arm[1,])] <- 1
    chr.arm.tab[2, i, which.max(chr.arm[2,])] <- 1
  }  
  
  return(chr.arm.tab)
}

AllelicGetCopyRatioComb <- function(q, delta, b, hscn.params) {
    xx <- (2 * delta * (c(1:q) - 1) + b)
    
    means <- Atten(xx, hscn.params[["fit.at"]])
    tx.means <- TxData(hscn.params, means)
    
    return(tx.means)
}

allelic_get_cr_grid_from_ccf_grid = function(qc, qs, delta, comb, ccf_grid) {
  d = qs - qc
  cr_grid = (2 * (delta * d) * ccf_grid) + comb[qc + 1]
  
  if (d < 0) {
    cr_grid = rev(cr_grid)
  }
  
  return(cr_grid)
}

GetAllelicChrArmSegs <- function(seg.obj, chr.arm.dat) {
    seg.dat <- seg.obj[["as.seg.dat"]]
    
    ##  chr, start, stop   
    uniq.seg.ix <- unique(seg.dat[, "seg.ix"])
    seg.ix <- seg.dat[, "seg.ix"]
    n.seg <- length(uniq.seg.ix)
    
    int.w <- c()
    low.ix <- c()
    high.ix <- c()
    
    arm.len.bp <- chr.arm.dat["End.bp"] - chr.arm.dat["Start.bp"]
    
    for (i in 1:n.seg) {
        if (sum(seg.ix == uniq.seg.ix[i]) != 2) {
            next
        }
        
        ix <- which(seg.ix == uniq.seg.ix[i])[1]
        
        if (as.integer(seg.dat[ix, "Chromosome"]) != as.integer(chr.arm.dat["chr"])) {
            next
        }
        
        int.start <- max(as.numeric(c(seg.dat[ix, "Start.bp"],
                                      chr.arm.dat["Start.bp"])))
        int.end <- min(as.numeric(c(seg.dat[ix, "End.bp"],
                                    chr.arm.dat["End.bp"])))
        
        if (int.start > int.end)  {
          next
        }  ## seg does not overlap region
        
        pair.ix <- which(seg.ix == uniq.seg.ix[i])
        
        int.len.bp <- int.end - int.start
        int.w <- c(int.w, int.len.bp / arm.len.bp)
        
        low.ix <- c(low.ix, pair.ix[which.min(seg.dat[pair.ix, "copy_num"])])
        high.ix <- c(high.ix, pair.ix[which.max(seg.dat[pair.ix, "copy_num"])])
    }
    
    int.w <- as.numeric(int.w)
    
    return(list(int.w=int.w, low.ix=low.ix, high.ix=high.ix))
}

AllelicGetAbsSegDat <- function(segobj, at) {
  Q <- 15
  qq <- Q
  
  ## Get column number of the max of each row and the expected 
  seg.qz.tab <- segobj[["mode.res"]][["seg.qz.tab"]][1, , ]
  seg.q.tab <- segobj[["mode.res"]][["seg.q.tab"]][1, , ]
  
  max.mat <- apply(seg.qz.tab, 1, which.max)
  subclonal.ix <- max.mat == (Q + 1)
  
  max.mat <- apply(seg.q.tab,1, which.max)
  exp.mat <- apply(seg.q.tab, 1, function(x) {
    x <- x[1:qq] / sum(x[1:qq])
    return(sum(x * c(1:qq)))
  })
  
  ## seg_list is relevant seg table
  seg.list <- segobj[["as.seg.dat"]]
  
  ## get unique seg_ix values
  u <- unique(seg.list[, "seg.ix"])
  
  ## make vectors of 0s for columns
  modal.a1 <- vector(mode="numeric", length=length(u))
  expected.a1 <- vector(mode="numeric", length=length(u))
  modal.a2 <- vector(mode="numeric", length=length(u))
  expected.a2 <- vector(mode="numeric", length=length(u))
  LOH <- vector(mode="numeric", length=length(u))
  HZ <- vector(mode="numeric", length=length(u))
  subclonal.a1 <- vector(mode="numeric", length=length(u))
  subclonal.a2 <- vector(mode="numeric", length=length(u))
  copy.ratio <- vector(mode="numeric", length=length(u))
  hscr.a1 <- vector(mode="numeric", length=length(u))
  hscr.a2 <- vector(mode="numeric", length=length(u))
  
  cancer.cell.frac.a1 <- vector(mode="numeric",length=length(u))
  cancer.cell.frac.a2 <- vector(mode="numeric",length=length(u))
  
  ccf.ci95.low.a1 <-  vector(mode="numeric",length=length(u))
  ccf.ci95.low.a2 <-  vector(mode="numeric",length=length(u))
  
  ccf.ci95.high.a1 <-  vector(mode="numeric",length=length(u))
  ccf.ci95.high.a2 <-  vector(mode="numeric",length=length(u))
  
  sc_tab = segobj$mode.res$subclonal_SCNA_res$subclonal_SCNA_tab[1, ,]
  ccf_hat = round(sc_tab[, "CCF_hat"], 5)
  ccf_ci95 = round(sc_tab[, c("CI95_low", "CI95_high")], 5)
  
  ## for each pair in seg_ix, print modal and expected values and check for HZ and LOH
  for (i in seq_along(u)) {
    seg <- u[i]
    usegs <- which(seg.list[, "seg.ix"] == seg)
    hscr <- InvAtten(sort(seg.list[usegs, "copy_num"]), at)
    copy.ratio[i] <- round(sum(hscr) / 2, 5)
    hscr.a1[i] <- round(hscr[1], 5)
    hscr.a2[i] <- round(hscr[2], 5)
    
    modal.a1[i] <- max.mat[usegs[1]] - 1
    modal.a2[i] <- max.mat[usegs[2]] - 1
    expected.a1[i] <- round(exp.mat[usegs[1]] - 1, 5)
    expected.a2[i] <- round(exp.mat[usegs[2]] - 1, 5)
    
    subclonal.a1[i] <- subclonal.ix[usegs[1]]
    subclonal.a2[i] <- subclonal.ix[usegs[2]]
    
    cancer.cell.frac.a1[i] = ccf_hat[usegs[1]]
    cancer.cell.frac.a2[i] = ccf_hat[usegs[2]]
    ccf.ci95.low.a1[i] = ccf_ci95[usegs[1], 1]
    ccf.ci95.high.a1[i] = ccf_ci95[usegs[1], 2]
    ccf.ci95.low.a2[i] = ccf_ci95[usegs[2], 1]
    ccf.ci95.high.a2[i] = ccf_ci95[usegs[2], 2]
    
    if ((modal.a1[i]==0) && (modal.a2[i]==0)) {
      HZ[i] <- 1
    } else if ((modal.a1[i] == 0) || (modal.a2[i]==0)) {
      LOH[i] <- 1
    }
  }       
  
  ## round and delete appropriate fields from existing seg table
  ix <- which(colnames(seg.list) %in% c("seg.ix", "copy.num"))
  tab <- round(seg.list[1:length(u), c(-ix)], 5)
  
  return(cbind(tab, copy.ratio, hscr.a1, hscr.a2, modal.a1, modal.a2,
               expected.a1, expected.a2, subclonal.a1, subclonal.a2,
               cancer.cell.frac.a1, ccf.ci95.low.a1, ccf.ci95.high.a1,
               cancer.cell.frac.a2, ccf.ci95.low.a2, ccf.ci95.high.a2,
               LOH, HZ))
}
