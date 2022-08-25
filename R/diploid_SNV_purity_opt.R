## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.


## The purpose of this function is to add additional purity/ploidy modes to the 
## initial set of provisional modes.   This optimization assumes genome_mass=2 
## and searches for the optimal alpha by fitting the somatic mutations.

run_diploid_snv_purity_opt = function(obs_scna, mut_cn_dat, pi_som_theta_q, 
                                      mut_class_w, dom, d_res = 0.1,
                                      verbose=FALSE) {
  seg = obs_scna$segtab
  ## expected copy-number, should be 1.0
  e_cr = sum(seg[, "W"] * seg[, "copy_num"])  
  
  ## remove muts that are not on diploid regs
  scna.ix = which(abs(seg[, "copy_num"] - e_cr) > 0.1)
  
  if (length(scna.ix) > 0) {
    ##  works for HSCR or total CR
    if ("A1.ix" %in% colnames(mut_cn_dat)) {
      nix = mut_cn_dat[, "A1.ix"] %in% scna.ix | mut_cn_dat[,"A2.ix"] %in% scna.ix
    } else {
      nix = mut_cn_dat[, "mut_seg_ix"] %in% scna.ix 
    }     
    
    mut_cn_dat = mut_cn_dat[!nix, , drop=FALSE]
  }
  
  if (nrow(mut_cn_dat) == 0) {
    if (verbose) {
      print("No muts left on unaltered genome")
    }
    return(NA)
  }
  
  comb_1d_ll = function(par, mut_cn_dat, dom) {
    alpha = par
    
    ## FIXME: This looks like it should be a ||
    if (alpha < dom[1] | alpha > dom[2]) {
      return( .Machine$double.xmax) 
    } 
    
    mode_info = alpha
    names(mode_info)="alpha"

    mult_res = dir_post_fit_somatic_multiplicity(mut_cn_dat, mode_info, mut_class_w, 
                                                 pi_som_theta_q, print_fit=FALSE)
    
    ll = sum(mult_res$post_Prs[, "LL"], na.rm=TRUE)  + dDirichlet(mult_res$som_theta_Q_MAP, 
                                                                  pi_som_theta_q, log.p=TRUE)
    return(-ll)
  }
  
  mut_cn_dat = cbind(mut_cn_dat, q_hat=2, HS_q_hat_1=1, HS_q_hat_2=1)
  
  d_grid = seq(dom[1], dom[2], length=1000)
  ll_grid = rep(NA, length(d_grid))
  for (i in seq_along(d_grid)) {
    ll_grid[i] = -comb_1d_ll(d_grid[i], mut_cn_dat, dom)
  }
  d_ll = diff(ll_grid) 
  sd_ll = d_ll > 0
  rs = rbind(sd_ll, c(0, sd_ll[1:(length(sd_ll - 1))]))
  zc = which(rs[1, ]== 0 & rs[2, ]==1)
  mode_ix = zc
  if (ll_grid[length(ll_grid)] > ll_grid[length(ll_grid) - 1]) { 
    mode_ix = c(mode_ix, length(ll_grid)) 
  }
  if (ll_grid[1] > ll_grid[2]) { 
    mode_ix = c(mode_ix, 1) 
  }
  mode_tab = array( NA, dim=c(length(mode_ix), 3))
  
  a = d_grid[mode_ix]
  tau = 2
  # solve so comb[1] == e_cr
  tau_p = (2 * (1 - a) + 2 * a) / (e_cr * a)  -  2 * (1 - a) / a
  
  ##  convert to tau corresponding to ploidy=2 for this sample
  b = 2 * (1 - a) / (2 * (1 - a) + a * tau_p)
  delta = a / (2 * (1 - a) + a * tau_p)
  
  mode_tab[, 1] = b
  mode_tab[, 2] = log(delta)
  mode_tab[, 3] = ll_grid[mode_ix]
     
  return(mode_tab)
}


