## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

total_get_muts_nearest_clonal_scna = function(mut_cn_dat, seg_qz_hat, seg_q_hat, Q) {
  if (! "mut_seg_ix" %in% colnames(mut_cn_dat)) { 
    stop("wrong colnames") 
  }
  
  muts_p_qz = which.max(seg_qz_hat[mut_cn_dat[, "mut_seg_ix"], ])
  sc_mat = which.max(muts_p_qz)  
  subclonal_ix = (sc_mat == Q)
  
  muts_p_q = array(NA, dim=c(nrow(mut_cn_dat), 1, Q))
  
  for (i in seq_len(nrow(mut_cn_dat))) {
    muts_p_q[i, 1, ] = seg_q_hat[mut_cn_dat[i, "mut_seg_ix"], ]   
  }
  
  # todo: integrate over this instead
  muts_q_hat = apply(muts_p_q[, 1, , drop=FALSE], 1, which.max ) - 1
  mut_cn_dat = cbind(mut_cn_dat,  q_hat=muts_q_hat, HS_q_hat_1=NA, HS_q_hat_2=NA) 
  
  return(mut_cn_dat) 
}

total_get_mut_seg_ix = function(maf, segtab) {
  ## compute lookup table for each mutation into seg_Q_tab
  N = nrow(maf)
  mut_seg_ix = matrix(NA, nrow=N, ncol=1)
  colnames(mut_seg_ix) = c("mut_seg_ix")
  
  ##  This is Dan-Avi's fault.
  if (!("End_position" %in% colnames(maf))) {
    end_position = maf[, "Start_position"]
    maf = cbind(maf, End_position=end_position)
  }
  
  for (i in seq_len(N)) {
    seg.ix = maf[i, "Chromosome"] == segtab[, "Chromosome"] &
      maf[i, "Start_position"] >= segtab[, "Start.bp"] &
      maf[i, "End_position"] <= segtab[, "End.bp"] 
    
    if (sum(seg.ix) != 1) { 
      next 
    }
    
    seg_id = which(seg.ix)
    mut_seg_ix[i, 1] = seg_id
  }
  
  # for debugging only.  In general, we want to propagate NAs for un-mappable muts and handle them correctly
  #   if( any(is.na(mut_seg_ix)) ) { stop() }
  
  return(mut_seg_ix)
}

get_tcn_somatic_mut_comb = function(alpha, q) {
  q_a = c(1:q)
  
  f_reads = alpha * q_a / (2 * (1 - alpha) + alpha * q)
  vals = q_a
  res = cbind(vals, f_reads)
  
  return(res)
}

get_tcn_germline_mut_comb = function(alpha, q) {
  eps = 1e-3
  q_vals = c(0:q)
  
  f_reads = ((1 - alpha) + alpha * q_vals) / (2 * (1 - alpha) + alpha * q)   
  f_reads[3] = 1 - eps
  vals = q_vals
  res = cbind(vals, f_reads)
  
  return(res)
}

# version for tCR SCNA data
## data are all from a single sample
total_calc_sample_muts_post_pr = function(mut_cn_dat, mode_info, mut_w, som_theta_q)  {
  cols = c("Pr_somatic_clonal", "Pr_germline", "Pr_subclonal", "Pr_subclonal_wt0", 
            "Pr_wt0", "Pr_ge2", "Pr_cryptic_SCNA", "LL")
  mut_pr_mat = matrix(NA, nrow=nrow(mut_cn_dat), ncol=length(cols))
  colnames(mut_pr_mat) = cols
  
  som_mut_q = matrix(NA, nrow=nrow(mut_cn_dat), ncol=length(som_theta_q))
  alpha = mode_info["alpha"]
    
  u_keys = unique(mut_cn_dat[, "q_hat"])
  n_keys = length(u_keys)
  mut_n_keys = rep(NA, n_keys)
  
  for (i in seq_len(n_keys)) {
    qt = u_keys[i]
    mut_ix = mut_cn_dat[, "q_hat"] == qt 
    mut_n_keys[i] = sum(mut_ix)
  }
  
  u_keys = u_keys[order(mut_n_keys, decreasing=TRUE)]
  mut_n_keys = sort(mut_n_keys, decreasing=TRUE)
  
  for (i in seq_len(n_keys)) {
    qt = u_keys[i]
    mut_ix = mut_cn_dat[, "q_hat"] == qt 
    
    ## muts on a HZ-del should reflect Pr of seq artifacts
    if (qt == 0) {
      mut_pr_mat[mut_ix, ] = 0
      mut_pr_mat[mut_ix, "LL"] = log(mut_w$OL)
      som_mut_q[mut_ix, ] = 0
      som_mut_q[mut_ix, 1] = 1
      next
    }
    
    ## Somatic model calc
    s_delta = alpha / (2 * (1 - alpha) + alpha * (qt))
    cov_vals = rowSums(mut_cn_dat[mut_ix, c("ref", "alt")])
    fhat_vals =  mut_cn_dat[mut_ix, "alt"] / cov_vals
    
    som_comb = get_tcn_somatic_mut_comb(alpha, qt)
    som_w= SomaticMutPrior( mut_w, nrow(som_comb), som_theta_q)
    som_w =  matrix(som_w, ncol=nrow(som_comb), nrow=length(fhat_vals), byrow=TRUE)
    som_log_pr = log(som_w) + FhatCombPost(fhat_vals, cov_vals, som_comb[, "f_reads"])
    
    som_pr = cbind(som_log_pr, matrix(-Inf, nrow=nrow(som_log_pr), 
                                      ncol=length(som_theta_q) - ncol(som_log_pr)))
    som_pr = som_pr - LogAdd(som_pr)
    som_mut_q[ mut_ix, ] = exp(som_pr)
    
    gl_w = GermlineMutPrior(mut_w)
    if (any(gl_w > 0)) {
      stop( "Germline multiplcity not implemented for tCR input" )
    } else { 
      germ_log_pr = matrix(log(0), ncol=length(gl_w), nrow=length(fhat_vals))
    }
    
    subclonal_log_pr = log(mut_w$SC) + SubclonalPost(fhat_vals, s_delta, cov_vals)
    cryptic_scna_log_pr = log(0.01) + CrypticScnaPost(fhat_vals, s_delta, cov_vals)
    
    mut_ll = cbind( germ_log_pr, cryptic_scna_log_pr, subclonal_log_pr, som_log_pr )
    z = LogAdd(mut_ll)
    pr = matrix(exp( mut_ll - z), nrow=nrow(mut_ll), ncol = ncol(mut_ll))
    mut_pr_mat[mut_ix, "LL"] = z
    
    if (any(!is.finite(z))) { 
      stop("Non-finite log likelihood") 
    }
    
    mut_pr_mat[mut_ix, "Pr_somatic_clonal"] = 1 - rowSums(pr[, 1:5, drop=FALSE])
    mut_pr_mat[mut_ix, "Pr_germline"] = rowSums(pr[, 1:3, drop=FALSE])
    mut_pr_mat[mut_ix, "Pr_subclonal"] = pr[, 5, drop=FALSE]
    mut_pr_mat[mut_ix, "Pr_cryptic_SCNA"] = pr[, 4, drop=FALSE]
    
    if (qt == 1) {
      mut_pr_mat[mut_ix, "Pr_wt0"] = pr[, ncol(pr), drop=FALSE]
      mut_pr_mat[mut_ix, "Pr_subclonal_wt0"] = mut_pr_mat[mut_ix, "Pr_subclonal"]
    }
    ## prob that a clonal somatic multiplicity is >= 2
    if  (qt > 1) {
      mut_pr_mat[mut_ix, "Pr_ge2"] = rowSums(pr[, 7:ncol(pr), drop=FALSE])
    }
  }
  
  ## assume subclonal muts have mult 1 in SC fraction.
  modal_q_s = rep(NA, nrow(som_mut_q))
  if  (mut_w$SM > 0) {
    nix =  mut_pr_mat[, "Pr_somatic_clonal"] < 0.5
    nix[is.na(nix)] = TRUE
    
    if (any(!nix)) {
      modal_q_s[!nix] = apply(som_mut_q[!nix, , drop=FALSE], 1, which.max)
    }
    modal_q_s[nix] = 1
  }
  
  post_prs = cbind(mut_pr_mat, modal_q_s=modal_q_s) 
  pr_somatic = rowSums(post_prs[, c("Pr_somatic_clonal", "Pr_subclonal"), drop=FALSE], na.rm=T)
  post_prs = cbind(post_prs, Pr_somatic=pr_somatic)
  
  return(list(post_Prs=post_prs, som_mut_Q_tab=som_mut_q))
}


mut_qm_fc_grid_integral = function(int_mat, outer_prod,  fhat, cov) {
  h_qm_ll_mat = matrix(NA, nrow=length(fhat), ncol=nrow(outer_prod))
  for (j in 1:length(fhat)) {
    h_qm_ll_mat[j, ] =  LogAdd(int_mat + dbeta(outer_prod, fhat[j] * cov[j] + 1, 
                                               (1 - fhat[j]) * cov[j] + 1, log=TRUE))
  }
  
  return(h_qm_ll_mat)
}

# fc is a vector 
delta_s = function(fc, qc, qs, alpha) {
  return(alpha / ( 2 * (1 - alpha) + alpha * qc * (1 - fc) + alpha * fc * qs))
}

# version for subclonal tCN SCNAs 
total_calc_sample_muts_on_subclonal_scna = function(mut_cn_dat, subclonal_scna_tab, 
                                                    scna_ccf_dens, mode_info, obs_scna, 
                                                    mut_w, som_theta_q) {
  ## data are all from a single sample
  cols = c( "Pr_somatic_clonal", "Pr_germline", "Pr_subclonal", "Pr_subclonal_wt0", "Pr_wt0", "Pr_ge2", "Pr_cryptic_SCNA", "LL" )
  mut_pr_mat = matrix( NA, nrow=nrow(mut_cn_dat), ncol=length(cols) )
  colnames(mut_pr_mat) = cols
  
  alpha = mode_info["alpha"]
  cov_vals = rowSums(mut_cn_dat[, c("ref", "alt")])
  fhat_vals =  mut_cn_dat[, "alt"] / cov_vals
  
  f_c_grid = as.numeric(colnames(scna_ccf_dens))
  d_f_c = c(0, diff(f_c_grid))
  f_c_hat = subclonal_scna_tab[, "CCF_hat"]
  
  for (i in seq_len(nrow(subclonal_scna_tab))) {
    if (!subclonal_scna_tab[i, "subclonal_ix"])  { 
      next 
    }
    
    mut_ix = mut_cn_dat[, "mut_seg_ix"] == i 
    n_mut = sum(mut_ix)
    
    if (n_mut == 0) { 
      next 
    }
    
    cat(paste(n_mut, ",", sep=""))
    
    fhat = fhat_vals[mut_ix]
    cov = cov_vals[mut_ix]
    
    f_c_dens = scna_ccf_dens[i, , drop=FALSE]
    qc = subclonal_scna_tab[i, "qc"] 
    qs = subclonal_scna_tab[i, "qs"] 
    
    ## can happen for segs < c0
    if (is.na(qc) || is.na(qs)) { 
      next 
    }   
    
    ## for all muts on seg i
    
    ## H1: subclonal SNV in disjoint subpop as subclonal SCNA
    qm_set = seq_len(qc)
    int_mat = matrix(log(f_c_dens) + log(d_f_c), nrow=length(qm_set), ncol=length(f_c_dens), byrow=TRUE)
    outer_prod = qm_set %*% t((1 - f_c_grid) * delta_s(f_c_grid, qc, qs, alpha))
    h1_qm_ll = mut_qm_fc_grid_integral(int_mat, outer_prod,  fhat, cov)
        
    # H2: subclonal SNV in same subpop as subclonal SCNA
    qm_set = seq_len(qs)
    int_mat = matrix(log(f_c_dens) + log(d_f_c), nrow=length(qm_set), ncol=length(f_c_dens), byrow=TRUE)
    outer_prod = qm_set %*% t( f_c_grid * delta_s(f_c_grid, qc, qs, alpha))
    h2_qm_ll = mut_qm_fc_grid_integral( int_mat, outer_prod,  fhat, cov )
    
    
    # H3: clonal SNV  
    # cis model
    qm_set = seq(max(1, qs - qc), qc)
    int_mat = matrix(log(f_c_dens) + log(d_f_c), nrow=length(qm_set), ncol=length(f_c_dens), byrow=TRUE)
    outer_prod = qm_set %*% t((1 - f_c_grid) * delta_s(f_c_grid, qc, qs, alpha)) + 
      (qm_set + (qs-qc))  %*%  t(f_c_grid * delta_s(f_c_grid, qc, qs, alpha))  
    h3_cis_qm_ll =  mut_qm_fc_grid_integral(int_mat, outer_prod,  fhat, cov)
    
    #  trans model
    if (qs > 0) {
      qm_set = seq(1, min(qc, qs))
      int_mat = matrix(log(f_c_dens) + log(d_f_c), nrow=length(qm_set), ncol=length(f_c_dens), byrow=TRUE)
      outer_prod = qm_set %*% t( delta_s(f_c_grid, qc, qs, alpha)) 
      h3_trans_qm_ll =  mut_qm_fc_grid_integral( int_mat, outer_prod,  fhat, cov )
    } else { 
      h3_trans_qm_ll = matrix(-Inf, nrow=n_mut, ncol=1)
    }
        
    ## unif SC and cryptic SCNA models
    sc_afs = c((1 - f_c_hat[i]) * delta_s(f_c_hat[i],  qc, qs, alpha), 
               f_c_hat[i] * delta_s(f_c_hat[i], qc, qs, alpha))
    max_sc_af = min(sc_afs)
    unif_sc_ll = rep(-Inf, n_mut)
    if (max_sc_af > 0) {
      for (j in seq_len(n_mut)) {
        unif_sc_ll[j] = UnifSubclonalPost(fhat[j], max_sc_af, cov[j]) 
      }
    }
    
    min_sc_af = max(sc_afs)
    cryptic_scna_ll = CrypticScnaPost(fhat, min_sc_af, cov)
    
    ## construct mixture weights, sum, and normalize likelyhoods
    qm_som_comb = get_tcn_somatic_mut_comb(alpha, max(qc, qs))
    
    norm_mut_q_w = function(h_qm_ll, som_theta_q, w) {
      n = ncol(h_qm_ll)
      res = som_theta_q[seq_len(n)]
      res = res / sum(res)
      qm_w = w * res
      
      # for all muts on seg
      som_w = matrix(qm_w, ncol=n, nrow=nrow(h_qm_ll), byrow=TRUE)
    }
    
    R = 1- (0.05 + mut_w$SC + 0.001)
    sc_scna = list(H1=0.05, H2= R/3, H3_cis = R/3, H3_trans = R/3,
                   unif_SC = mut_w$SC, cryptic_SCNA = 0.001)
        
    som_w = norm_mut_q_w(h1_qm_ll, som_theta_q, sc_scna$H1)
    h1_qm = log(som_w) + h1_qm_ll
    
    som_w = norm_mut_q_w(h2_qm_ll, som_theta_q, sc_scna$H2) 
    h2_qm = log(som_w) + h2_qm_ll 
    
    som_w = norm_mut_q_w(h3_cis_qm_ll, som_theta_q, sc_scna$H3_cis) 
    h3_cis_qm = log(som_w) + h3_cis_qm_ll
    
    som_w = norm_mut_q_w(h3_trans_qm_ll, som_theta_q, sc_scna$H3_trans) 
    h3_trans_qm = log(som_w) + h3_trans_qm_ll
    
    unif_SC = unif_sc_ll + log(mut_w$SC)
    cryptic_scna = cryptic_scna_ll + log(sc_scna$cryptic_SCNA)
    
    ##  complete likelihood for mixture model
    mut_ll = cbind(h1_qm, h2_qm, h3_cis_qm, h3_trans_qm, unif_SC, cryptic_scna)
    Z = LogAdd( mut_ll )
    if (any(!is.finite(Z))) { 
      stop("Non-finite log likelihood") 
    }
    
    h1_qm_pr = matrix(exp(h1_qm - Z), nrow=n_mut, ncol = ncol(h1_qm))
    h2_qm_pr = matrix(exp(h2_qm - Z), nrow=n_mut, ncol = ncol(h2_qm))
    h3_cis_qm_pr = matrix(exp(h3_cis_qm - Z), nrow=n_mut, ncol = ncol(h3_cis_qm))
    h3_trans_qm_pr = matrix(exp(h3_trans_qm - Z), nrow=n_mut, ncol = ncol(h3_trans_qm))
    pr_cryptic_scna = matrix(exp(cryptic_scna - Z), nrow=n_mut, ncol=1)
        
    ##  compute various marginals on posteriors for mut interpretation
    mut_pr_mat[mut_ix, "LL"] = Z
    pr_somatic_clonal = rowSums(cbind(h3_cis_qm_pr, h3_trans_qm_pr))
    mut_pr_mat[mut_ix, "Pr_somatic_clonal"] = pr_somatic_clonal
    # Not implemented yet
    mut_pr_mat[mut_ix, "pr_germline"] = NA  
    mut_pr_mat[mut_ix, "pr_subclonal"] =  1 - pr_somatic_clonal
    mut_pr_mat[mut_ix, "Pr_cryptic_SCNA"] = pr_cryptic_scna
    
    # wt0 == homozygous
    # Not implemented yet
    mut_pr_mat[mut_ix, "Pr_wt0"] = NA
    mut_pr_mat[mut_ix, "Pr_subclonal_wt0"] = NA
        
    ## prob that a clonal somatic multiplicity is > 1
    pr_ge2 = rep(0, n_mut)  
    if (ncol(h3_cis_qm_pr) > 1) {
      pr_ge2 = rowSums(h3_cis_qm_pr[, -1, drop=FALSE]) 
    }

    if (ncol(h3_trans_qm_pr) > 1) {
      pr_ge2 = pr_ge2 + rowSums(h3_trans_qm_pr[, -1, drop=FALSE]) 
    }
    mut_pr_mat[ mut_ix, "Pr_ge2"] = pr_ge2    
  }
  cat("\n")

  ##  Not implemented yet
  post_prs = cbind(mut_pr_mat, modal_q_s=NA)   
  
  pr_somatic = rowSums(post_prs[, c("Pr_somatic_clonal", "Pr_subclonal"), drop=FALSE], na.rm=T)
  post_prs = cbind(post_prs, Pr_somatic=pr_somatic)
    
  return(list(post_Prs=post_prs))
}


total_get_subclonal_scna_mut_ix = function(mut_cn_dat, subclonal_scna_tab) {
  return(mut_cn_dat[,"mut_seg_ix"] %in% which(subclonal_scna_tab[, "subclonal_ix"] == 1))
}
