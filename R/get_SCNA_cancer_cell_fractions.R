apply_subclonal_scna_model = function(segobj, mode_res, verbose=FALSE) {
  
  Q = dim(mode_res$theta.q.tab)[2]
  M = nrow(mode_res$mode.tab)
  n_seg = nrow(segobj$obs.scna$segtab)
  
  if (verbose) {
    print( paste("Evaluating subclonal SCNAs in ", M, " purity/ploidy modes: ", sep=""))
  }
  
  for (i in seq_len(M)) {
    res = get_scna_cancer_cell_fractions(segobj, mode_res, i)
    if (i == 1) {
      n_col = ncol(res$subclonal_SCNA_tab)
      subclonal_scna_tab = array(NA, dim=c(M, n_seg, n_col))
      dimnames(subclonal_scna_tab)[[3]] = colnames(res$subclonal_SCNA_tab)
            
      n_col = ncol(res$CCF_dens)
      ccf_dens = array(NA, dim=c(M, n_seg, n_col))
      dimnames(ccf_dens)[[3]] = colnames(res$CCF_dens)
    }
    
    subclonal_scna_tab[i, , ] = res$subclonal_SCNA_tab
    ccf_dens[i, , ] = res$CCF_dens
  }
  
  mode_res$subclonal_SCNA_res = list(subclonal_SCNA_tab = subclonal_scna_tab, CCF_dens = ccf_dens)
  
  return(mode_res)
}


get_scna_cancer_cell_fractions = function(segobj, mode_res, mode_ix, pr_subclonal_threshold=0.2) {
  mode_tab = mode_res$mode.tab 
  pr_subclonal = mode_res$seg.z.tab[mode_ix, ]
  seg_q_tab = mode_res$seg.q.tab[mode_ix, , ]
  seg_qz_tab = mode_res$seg.qz.tab[mode_ix, , ]
  
  obs = ExtractSampleObs(segobj)
  obs$error.model$fit.at = mode_tab[mode_ix, "AT"]
  n_seg = nrow(seg_q_tab)
  
  sigma_h = mode_tab[mode_ix, "sigma.h.hat"] 
  tau = mode_tab[mode_ix, "tau"]
  alpha = mode_tab[mode_ix, "alpha"]
  delta = alpha / (2 * (1 - alpha) + alpha * tau) 
  b = 2 * (1 - alpha) / (2 * (1 - alpha) + alpha * tau) 
  comb =  InvTxData(obs$error_model, GetCopyRatioComb(15, delta, b, obs$error_model))
  
  seg_qz_hat = apply(seg_qz_tab, 1, which.max) - 1
  modal_cn = which.max(colSums(seg_q_tab)) - 1
  scna_ix = seg_qz_hat != modal_cn
  
  ## convert dist over CR to dist over CCF for subclonal SCNAs
  subclonal_ix = pr_subclonal >= pr_subclonal_threshold
  res = get_subclonal_states(obs, subclonal_ix, modal_cn, comb)
  qc = res$qc
  qs = res$qs
  
  ccf_hat = rep(NA, n_seg)
  ccf_hat[!subclonal_ix] = 1
  
  ccf_ci95 = matrix(NA, nrow=n_seg, ncol=2)
  ccf_ci95[!subclonal_ix, ] = 1
  
  ccf_grid = seq( 0.01, 1, length=100)
  ccf_dens = matrix(NA, nrow=n_seg, ncol=length(ccf_grid))
  colnames(ccf_dens) = ccf_grid
  
  cr_dens = matrix(NA, nrow=n_seg, ncol=length(ccf_grid))
  cr_grid = matrix(NA, nrow=n_seg, ncol=length(ccf_grid))
  
  for (i in 1:nrow(ccf_dens)) {
    if (!subclonal_ix[i]) { 
      next 
      }
    if (is.na(qs[i])) { 
      next 
    }
    
    d = qs[i] - qc[i]
    
    ## 1st-order approx to attenuation effect
    dd = (comb[qs[i] + 1] - comb[qc[i] + 1])
    cr_grid[i, ] = (dd * ccf_grid) + comb[qc[i] + 1]
    
    cr_dens[i, ] = GetScnaStderrGridDensity(obs, cr_grid[i, ], sigma_h, i)    
    ccf_dens[i, ] = cr_dens[i, ] / sum(cr_dens[i, ])
    ccf_hat[i] = ccf_grid[which.max(ccf_dens[i, ])]
    ecdf = cumsum(ccf_dens[i, ])
    ccf_ci95[i, ] = approx(x=ecdf, y=ccf_grid, xout=c(0.025, 0.975))$y
  }
  
  nix1 = is.na(ccf_ci95[, 1])
  ccf_ci95[nix1, 1] = min(ccf_grid)
  nix2 = is.na(ccf_ci95[, 2])
  ccf_ci95[nix2, 2] = max(ccf_grid)
  
  ## Round up in last bin.   TODO round down in 1st bin 
  ix = ccf_ci95[, 2] > ccf_grid[length(ccf_grid) - 1]
  ccf_ci95[ix, 2] = 1.0
  colnames(ccf_ci95) = c("CI95_low", "CI95_high")
  
  subclonal_scna_tab = cbind(CCF_hat=ccf_hat, subclonal_ix=subclonal_ix, SCNA_ix=scna_ix, 
                             ccf_ci95, Pr_subclonal=pr_subclonal, qs=qs, qc=qc)
  res = list(subclonal_SCNA_tab=subclonal_scna_tab, CCF_dens=ccf_dens)
  
  return(res)
}


get_subclonal_states = function(obs, subclonal_ix, modal_cn, comb) {
  cr_set = InvTxData(obs$error_model, obs$d.tx)
  
  ## assume clonal CN is the modal value
  qc = rep(modal_cn, length(cr_set))  
  qs = rep(NA, length(cr_set))
  
  del_ix  = subclonal_ix & (cr_set < comb[modal_cn+1])
  gain_ix = subclonal_ix & (cr_set > comb[modal_cn+1])
  Q = length(comb)
  
  ## qs = CN in the subclone with the SCNA
  del_vals = (c(0:(modal_cn - 1)))
  for (i in seq_along(del_vals)) {
    qs[del_ix & (cr_set > comb[del_vals[i] + 1]) ] = del_vals[i]
  }
  
  ## dels below 0
  qs[del_ix & cr_set < comb[1]] = NA    
  
  qc[del_ix] = qs[del_ix] + 1
  
  gain_vals = rev(c((modal_cn + 1):Q))
  for (i in seq_along(gain_vals)) {
    qs[gain_ix & (cr_set < comb[gain_vals[i] + 1])] = gain_vals[i]
  }

  qc[gain_ix] = qs[gain_ix] - 1
  
  res = list(qs=qs, qc=qc)
  return(res)
}
