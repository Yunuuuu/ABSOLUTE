## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

FhatCombPost <- function(fhat, cov, f) {
  log.pr <- matrix(NA, nrow = length(fhat), ncol = length(f))
  for (i in 1:length(f)) {
    log.pr[, i] <- dbeta(f[i], fhat * cov + 1, (1 - fhat) * cov + 1, log = TRUE)
  }
  
  return(log.pr)
}

## FIXME: This is not used
ExpPowSubclonalPost <- function(f.hat.vals, sDelta, cov) {
  dd <- rep(NA, length(f.hat.vals))
  
  for (i in 1:length(f.hat.vals)) {
    MIN <- 1 / cov[i] 
    MAX <- 1 - MIN
    BY <- MIN
    grid <- seq(MIN, MAX, by=BY)

    exp.dens <- dexp(grid / sDelta, rate=10) * sDelta^-1
    exp.dens <- log(exp.dens / sum(exp.dens)) 
 
    s.dens <- power_lookup(n=cov[i], grid)
    s.dens <- log(s.dens / sum(s.dens)) + log(length(grid))

    dd[i] <- LogAdd(exp.dens + beta_dens + s.dens) +  log(length(grid))
  }

   return(dd)
}

ExpSubclonalPost <- function(f.hat.vals, sDelta, cov, lambda=25) {
  dd <- rep(NA, length(f.hat.vals))
  
  for (i in seq_along(f.hat.vals)) {
    MIN <- 1 / cov[i] 
    MAX <- 1 - MIN
    BY <- 0.01
    beta.grid <- seq(MIN, MAX, by=BY)
    
    exp.grid <- seq(MIN, MAX, by=BY)
    exp.dens <- dexp(exp.grid / sDelta, rate=lambda) * sDelta^-1
    exp.dens <- log(exp.dens / sum(exp.dens)) 
    
    beta.dens <- dbeta(beta.grid, f.hat.vals[i] * cov[i] + 1,
                       (1 - f.hat.vals[i]) * cov[i] + 1, log=TRUE)
    
    beta.dens <- beta.dens - LogAdd(beta.dens)
    
    dd[i] <- LogAdd(exp.dens + beta.dens) + log(length(beta.grid))
  }
  
  return(dd)
}

UnifSubclonalPost <- function(f.hat.vals, sDelta, cov, lambda=25) {
  max.sc.af <- sDelta
  
  loglik <- rep(NA, length(f.hat.vals))

  beta.int <- pbeta(max.sc.af, f.hat.vals * cov + 1, (1 - f.hat.vals) * cov + 1,
                    lower.tail=TRUE, log.p=TRUE)
  loglik <- beta.int + log(1 / max.sc.af)
  
  return(loglik)
}

CrypticScnaPost <- function(f.hat.vals, sDelta, cov) {
  min.sc.af <- sDelta

  loglik <- rep(NA, length(f.hat.vals))

  beta.int <- pbeta(min.sc.af, f.hat.vals * cov + 1, (1 - f.hat.vals) * cov + 1,
                    lower.tail=FALSE, log.p=TRUE)
  loglik <- beta.int + log(1 / min.sc.af)
  
  return(loglik)
}

SomaticMutPrior <- function(W, N, som.q.theta) {
  res <- som.q.theta[c(1:N)]
  res <- res / sum(res)
  return(W[["SM"]] * res)
}


GermlineMutPrior <- function(W) {
  res <- c(49, 49, 1)
  res <- res / sum(res)
  return(W[["GL"]] * res)
}

ClassifySomaticVariants <- function(prs, pr.thresh) {
  subclonal.ix <- prs[, "Pr_subclonal"] > pr.thresh
  subclonal.wt0.ix <- prs[, "Pr_subclonal_wt0"] > pr.thresh
  clonal.ix <- prs[, "Pr_somatic_clonal"] > pr.thresh
  wt0.ix <- prs[, "Pr_wt0"] > pr.thresh
  ge2.ix <- prs[, "Pr_ge2"] > pr.thresh

  clonal.het.ix <- prs[, "Pr_somatic_clonal"] * (1 - prs[, "Pr_wt0"]) > pr.thresh 
  homozygous.ix <- rowSums(prs[, c("Pr_wt0", "Pr_subclonal_wt0"), drop=FALSE]) > pr.thresh
    
  res <- cbind(subclonal.ix, subclonal.wt0.ix, clonal.ix, wt0.ix, clonal.het.ix, 
               ge2.ix, homozygous.ix)
  colnames(res) <- c("subclonal.ix", "subclonal_wt0.ix", "clonal.ix", "wt0.ix", 
                     "clonal_het.ix", "ge2.ix", "homozygous.ix")
  
  return(res)
}

## FIXME: Unused (I think this is future code)
GetGeneVarClassMat <- function(mut.dat, gene.list=NA, verbose=FALSE) {
  class.mat <- mut.dat[, c("subclonal.ix", "subclonal_wt0.ix", "clonal.ix",
                           "wt0.ix", "clonal_het.ix", "ge2.ix", "homozygous.ix")]
  
  ## remove unclassified mutations
  nix <- rowSums(class.mat[, c("subclonal.ix", "subclonal_wt0.ix",
                               "clonal.ix", "wt0.ix", "clonal_het.ix",
                               "ge2.ix")], na.rm=TRUE) == 0
  mut.dat <- mut.dat[!nix, ]
  class.mat <- class.mat[!nix, ]
  
  if (verbose) {
    print(paste(sum(nix), " unclassified variants removed", sep=""))
   }
  
  if (is.na(sum(nix))) {
    stop()
  }
  
  gene.names <- mut.dat[, "Hugo_Symbol"]
  
  ## Do all genes
  if (is.na(gene.list)) {
    gene.list <- unique(gene.names) 
  }
  
  cols <- c("Clonal het", "Clonal hom", "Clonal m >= 2", "Subclonal",
            "Subclonal hom", "N") 
  var.mat <- matrix(NA, ncol=length(cols), nrow=length(gene.list))
  rownames(var.mat) <- gene.list
  colnames(var.mat) <- cols
  
  for (i in 1:length(gene.list)) {
    gix <- (gene.names == gene.list[i])
    N <- sum(gix)
    
    var.mat[i, "Clonal het"] <- sum(gix & class.mat[, "clonal_het.ix"]) / N
    var.mat[i, "Subclonal"] <- sum(gix & class.mat[,  "subclonal.ix" ]) / N
    var.mat[i, "Clonal hom"] <- sum(gix & class.mat[, "wt0.ix"]) / N
    var.mat[i, "Clonal m >= 2"] <- sum(gix & class.mat[, "ge2.ix"]) / N
    var.mat[i, "Subclonal hom"] <- sum(gix & class.mat[, "subclonal_wt0.ix"]) / N
    var.mat[i, "N"] <- N
  }
  
  ix <- sort(rowSums(var.mat[, c("Clonal hom", "Clonal het")]),
             index.return=TRUE)$ix
  var.mat <- var.mat[ix, ]
  
  return(var.mat)
}

calc_ccf_posterior_grid = function(alt, ref, alpha, q, grid) {
  f = (alpha * grid) / (2 * (1 - alpha) + alpha * q)
  N = alt + ref
  
  ll_grid = dbinom(alt, N, f, log=TRUE )
  pr_grid = exp(ll_grid - LogAdd(ll_grid))
  
  return(pr_grid)   
}

calc_snv_ccf = function(mut_dat, mode_info, modal_q_s, pr_somatic_clonal) {
  muts = mut_dat
  alpha = mode_info["alpha"]
  cov = rowSums(muts[, c("alt", "ref")])
  f = muts[, "alt"] / cov
  Q = muts[, "q_hat"]
  
  q_s = rep(NA, length(pr_somatic_clonal))
  c.ix = pr_somatic_clonal > 0.5
  q_s[c.ix] = modal_q_s[c.ix]
  q_s[!c.ix] = 1 
  
  som_delta = alpha / (2 * (1 - alpha) + alpha * Q)
  
  f_s = q_s * som_delta    
  mut_scale = 1 / f_s   
  cell_frac = f * mut_scale
  cell_mult = f / som_delta
  
  ccf_grid = seq(0, 1, by=0.01)
  ccf_dens = matrix(NA, nrow=nrow(mut_dat), ncol=length(ccf_grid))
  ccf_ci95 = matrix(NA, nrow=nrow(ccf_dens), ncol=2)
  ccf_hat = rep(NA, nrow(mut_dat))
  
  for (i in seq_len(nrow(mut_dat))) {
    if (mut_dat[i, "q_hat"] == 0) { 
      next 
    }
    
    ccf_dens[i, ] = calc_ccf_posterior_grid(mut_dat[i, "alt"], mut_dat[i, "ref"], 
                                            alpha, mut_dat[i, "q_hat"], ccf_grid)
    ccf_hat[i] = ccf_grid[which.max(ccf_dens[i, ])]
    
    ecdf = cumsum(ccf_dens[i, ])

    ccf_ci95[i, ] = approx(x=ecdf, y=ccf_grid, xout=c(0.025, 0.975))$y
  }
  nix1 = is.na(ccf_ci95[, 1])
  ccf_ci95[nix1, 1] = min(ccf_grid)
  nix2 = is.na(ccf_ci95[,2])
  ccf_ci95[nix2, 2] = max(ccf_grid)
  
  ## Round up in last bin.   TODO round down in 1st bin 
  ix = ccf_ci95[, 2] > ccf_grid[length(ccf_grid) - 1]
  ccf_ci95[ix, 2] = 1.0
  
  ix = mut_dat[i, "q_hat"] == 0
  ccf_ci95[ix, ] = NA
  
  res = cbind(cell_mult, cell_frac, ccf_hat, ccf_ci95)
  colnames(res) = c("cell_mult", "old_cancer_cell_frac", "cancer_cell_frac", 
                    "ccf_CI95_low", "ccf_CI95_high")
  
  return(res)
}


dir_post_fit_somatic_multiplicity = function(mut_cn_dat, mode_info, mut_class_w, 
                                             pi_som_theta_q, print_fit=TRUE) {
  ## Dir mode
  som_theta_Q_mode = (pi_som_theta_q - 1) / (sum(pi_som_theta_q) - length(pi_som_theta_q)) 
  
  ## mut_class_w Dir mode
  mcv = c(mut_class_w$Pi_SM, mut_class_w$Pi_SC)
  mcw = (mcv - 1) / (sum(mcv) - length(mcv))   
  mut_class_w$SM = mcw[1]
  mut_class_w$SC = mcw[2]
  
  mult_res = CalcSampleMutsPostPr(mut_cn_dat, mode_info, mut_class_w, som_theta_Q_mode)
  
  mut_mat =  mult_res$som_mut_Q_tab * matrix(mult_res$post_Prs[,"Pr_somatic_clonal"], 
                                             nrow=nrow(mult_res$som_mut_Q_tab), 
                                             ncol=ncol(mult_res$som_mut_Q_tab), 
                                             byrow=TRUE) 
  
  nix = mult_res$post_Prs[, "Pr_somatic_clonal"] == 0
  mut_mat = mut_mat[!nix, , drop=FALSE]
  if (nrow(mut_mat) > 0) {
    som_q = colSums(mut_mat, na.rm=TRUE)
  } else { 
    som_q = rep(0, length(pi_som_theta_q)) 
  }
  
  ## Dir mode
  som_theta_q_map = (pi_som_theta_q + som_q - 1) / (sum(pi_som_theta_q + som_q) - length(pi_som_theta_q)) 
  ## FIXME: verbose
  if (print_fit) {  
    cat("som_theta_Q_MAP: ")
    print(round(som_theta_q_map, 5)) 
  }
  
  ## learn mixture weights and re-run
  mut_w =  mult_res$post_Prs[, c("Pr_somatic_clonal", "Pr_subclonal"), drop=FALSE]
  mut_w = mut_w / rowSums(mut_w)
  mut_pr = colSums(mut_w, na.rm=TRUE)
  
  ## mut_class_w Dir mode
  mcw = (mcv + mut_pr - 1) / (sum(mcv + mut_pr) - length(mcv))   
  mut_class_w$SM = mcw[1]
  mut_class_w$SC = mcw[2]
  
  if (print_fit) { 
    print(mut_class_w) 
  }
  
  ## Clonal muts only
  ## data are all from a single sample
  mult_res = CalcSampleMutsPostPr(mut_cn_dat, mode_info, mut_class_w, som_theta_q_map)   
  mult_res$som_theta_Q_MAP = som_theta_q_map
  
  return(mult_res) 
}

