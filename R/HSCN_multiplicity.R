## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

GetHscnSomaticMutComb <- function(alpha, q1, q2) {
  q.a <- c(1:q2)
  q <- q1 + q2
  f.reads <- alpha * q.a/(2 * (1 - alpha) + alpha * q)
  vals <- q.a
  res <- cbind(vals, f.reads)
  return(res)
}

GetHscnGermlineMutComb <- function(alpha, q1, q2) {
  eps <- 0.001
  
  q.vals <- c(q1, q2, q1 + q2)
  Q <- q1 + q2
  f.reads <- ((1 - alpha) + alpha * q.vals)/(2 * (1 - alpha) + alpha * Q)
  f.reads[3] <- 1 - eps
  vals <- q.vals
  res <- cbind(vals, f.reads)
  
  return(res)
}

AllelicGetMutSegIx <- function(maf, segtab) {
  ## compute lookup table for each mutation into seg.q.tab
  mut.seg.ix <- matrix(NA, nrow=nrow(maf), ncol=2)
  colnames(mut.seg.ix) <- c("A1.ix", "A2.ix")

  ## This is Dan-Avi's fault
  if (!("End_position" %in% colnames(maf))) {
    end.position <- maf[, "Start_position"]
    maf <- cbind(maf, "End_position"=end.position)
  }
  
  for (i in seq_len(nrow(maf))) {
    seg.ix <- (maf[ i, "Chromosome" ] == segtab[,"Chromosome"]) &
              (maf[ i, "Start_position"] >= segtab[,"Start.bp"]) &
              (maf[ i, "End_position"] <= segtab[,"End.bp"])

    if (sum(seg.ix) != 2) {
      next
    }
    
    seg.ids <- which(seg.ix)
    mut.seg.ix[i, 1] <- seg.ids[1] 
    mut.seg.ix[i, 2] <- seg.ids[2] 
  }
  
  ## for debugging only.  In general, we want to propagate NAs for un-mappable muts
  ## and handle them correctly
  return(mut.seg.ix)
}

allelic_get_muts_nearest_clonal_scna <- function(mut.cn.dat, seg.qz.hat, seg.q.hat, Q) {
  if (!all(c("A1.ix", "A2.ix") %in% colnames(mut.cn.dat)) ) {
    stop("wrong colnames")
  }
  
  muts.p.qz <- matrix(NA, nrow=nrow(mut.cn.dat), ncol=2)
  muts.p.qz[, 1] <- which.max(seg.qz.hat[ mut.cn.dat[, "A1.ix"], ])
  muts.p.qz[, 2] <- which.max(seg.qz.hat[ mut.cn.dat[, "A2.ix"], ]) 

  sc.mat <- cbind(which.max(muts.p.qz[, 1]),  which.max(muts.p.qz[, 2]))
  subclonal.ix <- apply((sc.mat == Q), 1, any)
  
  ## remove muts at subclonal SCNA regions for likelihood calc
  mut.cn.dat <- mut.cn.dat[!subclonal.ix, ]
  muts.p.q <- array(NA, dim=c( nrow(mut.cn.dat), 2, Q))
  
  for (i in seq_len(nrow(mut.cn.dat))) {
    muts.p.q[i, 1, ] <- seg.q.hat[mut.cn.dat[i,"A1.ix"], ]   
    muts.p.q[i, 2, ] <- seg.q.hat[mut.cn.dat[i,"A2.ix"], ]   
  }
  
  ## todo: integrate over this instead
  muts.q.hat <- cbind(apply(muts.p.q[, 1, ], 1, which.max ) - 1,
                      apply(muts.p.q[, 2, ], 1, which.max ) - 1)

  ## needed for following calc
  mut.cn.dat <- cbind(mut.cn.dat,  "q_hat" = rowSums(muts.q.hat),
                      "HS_q_hat_1"=muts.q.hat[,1],
                      "HS_q_hat_2"=muts.q.hat[,2] ) 
  
  return(mut.cn.dat) 
}

AllelicCalcSampleMutsPostPr <- function(mut.cn.dat, mode_info, mut.w, som.theta.q) {
  ## data are all from a single sample
  cols <- c("Pr_somatic_clonal", "Pr_germline", "Pr_subclonal",
            "Pr_subclonal_wt0", "Pr_wt0", "Pr_ge2", "Pr_GL_som_HZ_alt",
            "Pr_GL_som_HZ_ref", "Pr_cryptic_SCNA", "LL")
  mut.pr.mat <- matrix(NA, nrow=nrow(mut.cn.dat), ncol=length(cols))
  colnames(mut.pr.mat) <- cols
  
  som.mut.q <- matrix(NA, nrow=nrow(mut.cn.dat), ncol=length(som.theta.q))
  ## alpha <- mut.cn.dat[1, "purity"]
  alpha = mode_info["alpha"]
  hap.q.keys <- paste(mut.cn.dat[, "HS_q_hat_1"], mut.cn.dat[, "HS_q_hat_2"],
                      sep="__" )
  
  u.keys<- unique(hap.q.keys)
  n.keys <- length(u.keys)
  mut.n.keys <- rep(NA, n.keys)
  
  for(i in seq_len(n.keys)) {
    res <- strsplit(u.keys[i], "__" )[[1]]
    q1 <- as.integer(res[1])
    q2 <- as.integer(res[2])
    mut.ix <- (mut.cn.dat[, "HS_q_hat_1"] == q1) & (mut.cn.dat[, "HS_q_hat_2"] == q2)

    mut.n.keys[i] <- sum(mut.ix)
  }
  
  u.keys <- u.keys[ order(mut.n.keys, decreasing=TRUE) ]
  mut.n.keys <- sort(mut.n.keys, decreasing=TRUE)
  
  for (i in seq_len(n.keys)) {
    res <- strsplit(u.keys[i], "__" )[[1]]
    q1 <- as.integer(res[1])
    q2 <- as.integer(res[2])
    mut.ix <- (mut.cn.dat[, "HS_q_hat_1"] == q1) & (mut.cn.dat[, "HS_q_hat_2"] == q2) 

    ## FIXME: Is this really a &&?
    if ((q1 == 0) & (q2 == 0)) {
      mut.pr.mat[mut.ix, ] <- 0
      mut.pr.mat[mut.ix, "LL"] <- log(mut.w[["OL"]])
      som.mut.q[mut.ix, ] <- 0
      som.mut.q[mut.ix, 1] <- 1
      next
    }

    ## Somatic model calc
    sDelta <- alpha / (2 * (1 - alpha) + alpha * (q1 + q2))
    cov.vals <- rowSums(mut.cn.dat[mut.ix, c("ref", "alt")])
    fhat.vals <-  mut.cn.dat[mut.ix, "alt"] / cov.vals
    
    som.comb <- GetHscnSomaticMutComb(alpha, q1, q2)
    som.w <- SomaticMutPrior(mut.w, nrow(som.comb), som.theta.q)
    som.w <-  matrix(som.w, ncol=nrow(som.comb), nrow=length(fhat.vals), byrow=TRUE)
    som.log.pr <- log(som.w) + FhatCombPost(fhat.vals,
                                            cov.vals, som.comb[,"f.reads"] )
    som.pr = cbind(som.log.pr, matrix(-Inf, nrow=nrow(som.log.pr), 
                                      ncol=length(som.theta.q) - ncol(som.log.pr)))
    som.pr = som.pr - LogAdd(som.pr)
    som.mut.q[mut.ix, ] <- exp(som.pr)

    ## Germline model calc
    gl.w <- GermlineMutPrior(mut.w)
    if (any(gl.w > 0)) {
      germ.comb <- GetHscnGermlineMutComb(alpha, q1, q2)
      gl.w <-  matrix(gl.w, ncol=nrow(germ.comb), nrow=length(fhat.vals), byrow=TRUE) 
      germ.log.pr <- log(gl.w) +
        FhatCombPost(fhat.vals, cov.vals, germ.comb[,"f.reads"])
    } else { 
      germ.log.pr <- matrix(log(0), ncol=length(gl.w), nrow=length(fhat.vals))
    }

    ## Subclonal model calc
    subclonal.log.pr <- log(mut.w[["SC"]]) +
      SubclonalPost(fhat.vals, sDelta, cov.vals)
    cryptic.scna.log.pr <- log(0.01) + CrypticScnaPost(fhat.vals, sDelta, cov.vals)
    mut.ll <- cbind(germ.log.pr, cryptic.scna.log.pr, subclonal.log.pr, som.log.pr)
    Z <- LogAdd(mut.ll)
    Pr <- matrix(exp(mut.ll - Z), nrow=nrow(mut.ll), ncol = ncol(mut.ll))
    mut.pr.mat[mut.ix, "LL"] <- Z
    
    if (any(!is.finite(Z))) {
      stop("Non-finite log likelihood")
    }
    
    mut.pr.mat[mut.ix, "Pr_somatic_clonal"] <- 1 - rowSums(Pr[, c(1:5), drop=FALSE])
    mut.pr.mat[mut.ix, "Pr_germline"] <- rowSums(Pr[, c(1:3), drop=FALSE])
    mut.pr.mat[mut.ix, "Pr_subclonal"] <- Pr[, 5, drop=FALSE]
    mut.pr.mat[mut.ix, "Pr_cryptic_SCNA"] <- Pr[, 4, drop=FALSE]
    
    ## FIXME: This look like it should be a &&
    if ((q1 == 0) & (q2 > 0)) {
      mut.pr.mat[mut.ix, "Pr_wt0"] <- Pr[, ncol(Pr), drop=FALSE]
      ## probability the event is germline and somatically homozygous
      mut.pr.mat[mut.ix, "Pr_GL_som_HZ_alt"] <- Pr[, 2] / rowSums(Pr[, 1:3, drop=FALSE])
      mut.pr.mat[mut.ix, "Pr_GL_som_HZ_ref"] <- Pr[, 1] / rowSums(Pr[, 1:3, drop=FALSE])
    } else{
      mut.pr.mat[mut.ix, "Pr_wt0"] <- 0
    }
    
    if ((q1 == 0) & (q2 == 1)) {
      mut.pr.mat[mut.ix, "Pr_subclonal_wt0"] <- mut.pr.mat[mut.ix, "Pr_subclonal"]
    } else{
      mut.pr.mat[mut.ix, "Pr_subclonal_wt0"] <- 0
    }
    
    ## prob that a clonal somatic multiplicity is >= 2
    if (any(c(q1,q2) > 1)) {
      mut.pr.mat[mut.ix, "Pr_ge2"] <- rowSums(Pr[, c(7 : ncol(Pr)), drop=FALSE])
    }
  }
  
  ## assume subclonal muts have mult 1 in SC fraction.
  modal.q.s <- rep(NA, nrow(som.mut.q))
  if (mut.w[["SM"]] > 0) {
    nix <- mut.pr.mat[, "Pr_somatic_clonal"] < 0.5
    nix[is.na(nix)] = TRUE
    if (any(!nix)) {
      modal.q.s[!nix] <- apply(som.mut.q[!nix, , drop=FALSE], 1, which.max)
    }
    modal.q.s[nix] <- 1
  }
  
  post.prs <- cbind(mut.pr.mat, modal_q_s=modal.q.s)
  pr.somatic <- rowSums(post.prs[, c("Pr_somatic_clonal", "Pr_subclonal"), drop=FALSE], na.rm=TRUE)
  post.prs <- cbind(post.prs, Pr_somatic=pr.somatic)
  
  return(list(post_Prs=post.prs, som_mut_Q_tab=som.mut.q))
}

## data are all from a single sample
allelic_calc_sample_muts_on_subclonal_scna = function(mut_cn_dat, mode_info, obs_scna) {
  
  cols = c("Pr_somatic_clonal", "Pr_germline", "Pr_subclonal", "Pr_subclonal_wt0", 
           "Pr_wt0", "Pr_ge2", "Pr_GL_som_HZ_alt", "Pr_GL_som_HZ_ref", "Pr_cryptic_SCNA", "LL")
  mut_pr_mat = matrix(NA, nrow=nrow(mut_cn_dat), ncol=length(cols))
  colnames(mut_pr_mat) = cols
  
  alpha = mode_info["alpha"]
  cov_vals = rowSums( mut_cn_dat[, c("ref", "alt")])
  fhat_vals =  mut_cn_dat[, "alt"] / cov_vals
    
  # loop over HSCR segs
  a_ix_keys = paste(mut_cn_dat[, "A1.ix"], mut_cn_dat[, "A2.ix"], sep="__" )
  u_keys= unique(a_ix_keys)
  n_keys = length(u_keys)
  mut_n_keys = rep(NA, n_keys)
  
  for (i in 1:n_keys) {
    res = strsplit(u_keys[i], "__" )[[1]]
    A1.ix = as.integer(res[1])
    A2.ix = as.integer(res[2])
    mut_ix = mut_cn_dat[, "A1.ix"] == A1.ix & mut_cn_dat[,"A2.ix"] == A2.ix 
    mut_n_keys[i] = sum(mut_ix)
  }
  
  u_keys = u_keys[order(mut_n_keys, decreasing=TRUE)]
  mut_n_keys = sort(mut_n_keys, decreasing=TRUE)
  
  cr_grid = seq( 0, 5, length=1000 )
  cr_dens = get_scna_stderr_grid_density( obs_scna, cr_grid, mode_info["sigma_H_hat"] )
  cr_dens = cr_dens / rowSums(cr_dens)
  
  cn_grid = seq(0,5, length=1000)
  scale = 1 / mode_info["Delta"]
  shift = -mode_info["b"] * scale
  x = cr_grid * scale + shift
  
  for(i in 1:nrow(cr_dens)) {
    y = cr_dens[i, ] * abs(scale)^-1 
    cn_dens[i, ] = approx( x, y, xout=cn_grid)$y
  }      
  
  stop()
  ## FIXME: WTH?
  
  for( i in 1:n_keys )
  {
    res = strsplit(u_keys[i], "__" )[[1]]
    A1.ix = as.integer(res[1])
    A2.ix = as.integer(res[2])
    mut_ix = mut_cn_dat[,"A1.ix"] == A1.ix & mut_cn_dat[,"A2.ix"] == A2.ix 
    
        
    pr_f1 =  (cn_dens[A1.ix,] * tcn_dens ) 
    
    
    f1 =  (cn_grid * alpha) / ( 2*(1-alpha) + alpha*cn_grid ) 
    
    H1_LL = sum(  pr_f * 
      dbeta( f1, fhat*cov + 1, (1-fhat)*cov + 1, log=TRUE ) )
    
  }
}


allelic_get_subclonal_scna_mut_ix = function(mut_cn_dat, subclonal_scna_tab) {
  ix = mut_cn_dat[, "A1.ix"] %in% which(subclonal_scna_tab[, "subclonal_ix"] == 1) | 
    mut_cn_dat[, "A2.ix"] %in% which(subclonal_scna_tab[, "subclonal_ix"] == 1)
  
  return(ix)
}
