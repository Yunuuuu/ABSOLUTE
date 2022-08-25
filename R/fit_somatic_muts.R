## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

ApplySomaticMutsModel <- function(mode.res, obs_scna, mut.cn.dat, pi.som.theta.q,
                                  mut_class_w, Q, verbose=FALSE) {
  q <- dim(mode.res[["theta.q.tab"]])[2]

  if (verbose) {
    print(paste("Evaluating ", nrow(mut.cn.dat),
                " mutations over ", nrow(mode.res[["mode.tab"]]),
                " purity/ploidy modes: ", 
                sep = ""))
  }
  
  for (j in 1:nrow(mode.res[["mode.tab"]])) {
    alpha <- mode.res[["mode.tab"]][j, "alpha"]
    seg.qz.hat <- mode.res[["seg.qz.tab"]][j, , ]
    seg.q.hat <- mode.res[["seg.q.tab"]][j, , ]
    
    ##  sample purity     
    subclonal_scna_tab = mode.res$subclonal_SCNA_res$subclonal_SCNA_tab[j, , ]
    ccf_dens = mode.res$subclonal_SCNA_res$CCF_dens[j, , ]
    
    res = FitSomaticMuts(mut.cn.dat, mode.res$mode.tab[j, ], subclonal_scna_tab, 
                         ccf_dens, obs_scna, seg.qz.hat, seg.q.hat, pi.som.theta.q, 
                         mut_class_w, Q, verbose=verbose)
        
    muts.post.prs <- res[["muts.post.prs"]]
    som.theta.q.map <- res[["som.theta.q.map"]]
    
    if (j == 1) {
      mode.res[["muts.post.prs"]] <- array(NA, dim = c(nrow(mut.cn.dat),
                                                 ncol(muts.post.prs) + 1,
                                                 nrow(mode.res[["mode.tab"]])))
      dimnames(mode.res[["muts.post.prs"]])[[2]] <- c(colnames(muts.post.prs), "purity")
      dimnames(mode.res[["muts.post.prs"]])[[1]] <- rownames(muts.post.prs)
    }
    
    mode.res[["muts.post.prs"]][, , j] <- cbind(muts.post.prs, purity = alpha)
    
    mode.res[["mode.tab"]][j, "somatic.mut.ll"] <- sum(muts.post.prs[, "LL"],
                                                       na.rm = TRUE) + 
                                                         dDirichlet(som.theta.q.map,
                                                                    pi.som.theta.q,
                                                                    log.p = TRUE)
    if (verbose) {
      cat(".")
    }
  }
  if (verbose) {
    cat("\n")
  }
  
  new.ll <- mode.res[["mode.tab"]][, "clust_LL"] +
    mode.res[["mode.tab"]][, "log.partition"] + 
      mode.res[["mode.tab"]][, "somatic.mut.ll"]
  mode.res[["mode.tab"]][, "post_LL"] <- new.ll
  
  return(mode.res)
}

FitSomaticMuts <- function(mut.cn.dat, mode_info, subclonal_scna_tab, scna_ccf_dens,
                           obs_scna, seg.qz.hat, seg.q.hat, pi.som.theta.q, 
                           mut_class_w, Q, verbose=FALSE) {
  
  mut.cn.dat = get_muts_nearest_clonal_scna(mut.cn.dat, seg.qz.hat, seg.q.hat, Q)
  
  mult_res = dir_post_fit_somatic_multiplicity(mut.cn.dat, mode_info, mut_class_w, pi.som.theta.q,
                                               print_fit=verbose)
  som_theta_q_map = mult_res$som_theta_Q_MAP
  post_prs = mult_res$post_Prs
  
  clonal_scna_mut_ix = !get_subclonal_scna_mut_ix(mut.cn.dat, subclonal_scna_tab)
  
  ## FIXME: next version - not fully implemented yet
  if (FALSE)  {
    # For muts on clonal SCNAs
    clonal_scna_mult_res = dir_post_fit_somatic_multiplicity(mut.cn.dat[clonal_scna_mut_ix, ], 
                                                             mode_info, mut_class_w, pi.som.theta.q,
                                                             print_fit=verbose)
    som_theta_q_map = clonal_scna_mult_res$som_theta_Q_MAP
    post_prs = matrix(NA, nrow=nrow(mut.cn.dat), ncol=ncol(clonal_scna_mult_res$post_prs))
    colnames(post_prs) = colnames(clonal_scna_mult_res$post_Prs)
    post_prs[clonal_scna_mut_ix, ] = clonal_scna_mult_res$post_Prs 
    
    # " on subclonal SCNAs
    subclonal_scna_mult_res = calc_sample_muts_on_subclonal_scna(mut.cn.dat[!clonal_scna_mut_ix, ], 
                                                                 subclonal_scna_tab, 
                                                                 scna_ccf_dens, mode_info, obs_scna, 
                                                                 mut_class_w, som_theta_q_map)
    post_prs[!clonal_scna_mut_ix, ] = subclonal_scna_mult_res$post_Prs 
  }
    
  ## Subclonal SCNA?
  var_classes = ClassifySomaticVariants(post_prs, 0.5)
  q_s = post_prs[, "modal_q_s"]
  mut_ccf_res = calc_snv_ccf(mut.cn.dat, mode_info, q_s, post_prs[, "Pr_somatic_clonal"])
  ##
  
  muts.post.prs <- cbind(as.matrix(mut.cn.dat[, c("q_hat", "HS_q_hat_1",
                                                  "HS_q_hat_2"), drop=FALSE]),
                         post_prs, var_classes, mut_ccf_res, 
                         subclonal_SCNA=!clonal_scna_mut_ix)
                                                    
  return(list(muts.post.prs = muts.post.prs, som.theta.q.map=som_theta_q_map))
}

CreateMutCnDat <- function(maf, seg.dat, min.mut.af, verbose=FALSE) {
  ## todo - add this for preprocessing in Multiplicity/
  mut.cn.dat <- maf
  
  if ("total_normals_called" %in% colnames(mut.cn.dat)) {
    ix <- mut.cn.dat[, "total_normals_called"] > 1
    if (verbose) {
      print(paste("Removing ", sum(ix), " of ", length(ix),
                  " mutations due to seen in > 1 normals", 
                  sep = ""))
    }
    mut.cn.dat <- mut.cn.dat[!ix, ]
  }
  
  if ("dbSNP_Val_Status" %in% colnames(mut.cn.dat)) {
    mut.cn.dat[["dbSNP_Val_Status"]][is.na(mut.cn.dat[["dbSNP_Val_Status"]])] <- ""
  }
  
  cols <- colnames(mut.cn.dat)
  
  cix <- which(cols %in% c("i_t_ref_count", "t_ref_count"))
  if (length(cix) > 0) {
    colnames(mut.cn.dat)[cix] <- "ref"
  } else {
    stop("Malformed MAF file, no ref column supplied")
  }

  cix <- which(cols %in% c("i_t_alt_count", "t_alt_count"))
  if (length(cix) > 0) {
    colnames(mut.cn.dat)[cix] <- "alt"
  } else {
    stop("Malformed MAF file, no alt column supplied")
  }
  
  cix <- which(cols == "dbSNP_Val_Status")
  if (length(cix) > 0) {
    colnames(mut.cn.dat)[cix] <- "dbSNP"
  } else {
    stop("Malformed MAF file, no dbSNP_Val_Status column supplied")
  }
  
  cix <- which(cols == "Tumor_Sample_Barcode")
  if (length(cix) > 0) {
    colnames(mut.cn.dat)[cix] <- "sample"
  } else {
    stop("Malformed MAF file, no Tumor_Sample_Barcode column supplied")
  }
  
  na.ix <- apply(is.na(mut.cn.dat[, c("ref", "alt")]), 1, sum) > 0
  if (verbose) {
    print(paste("Removing ", sum(na.ix), " of ", length(na.ix),
                " mutations with NA coverage", 
                sep = ""))
  }
  mut.cn.dat <- mut.cn.dat[!na.ix, ]
  
  af <- mut.cn.dat[, "alt"] / (mut.cn.dat[, "alt"] + mut.cn.dat[, "ref"])
  ix <- af < min.mut.af
  if (verbose) {
    print(paste("Removing ", sum(ix), " of ", length(ix),
                " mutations due to allelic fraction < ", 
                min.mut.af, sep = ""))
  }
  mut.cn.dat <- mut.cn.dat[!ix, ]
  
  # remove sites in IG regions
  ix = mut.cn.dat[, "Hugo_Symbol"] %in% c("ADAM6")
  if (verbose) {
    print(paste("Removing ", sum(ix), " mutations in IG regions", sep=""))
  }
  if (sum(!ix) == 0) {  
    stop("no mutations left!") 
  }
  mut.cn.dat = mut.cn.dat[!ix, , drop=FALSE]
    
  mut.seg.ix <- GetMutSegIx(mut.cn.dat, seg.dat[["obs.scna"]][["segtab"]])  
  ix <- apply(is.na(mut.seg.ix), 1, sum) == 0
  
  if (verbose && (sum(!ix) > 0)) {
     
    print(paste("Removing ", sum(!ix), " unmapped mutations on Chrs: ", sep = ""))
    print(mut.cn.dat[!ix, "Chromosome"])
  }
  if (sum(ix) == 0) {
    stop("No mutations left")
  }
  mut.cn.dat <- mut.cn.dat[ix, ]
  mut.seg.ix <- mut.seg.ix[ix, , drop = FALSE]
  mut.cn.dat <- cbind(mut.cn.dat, mut.seg.ix)
  return(mut.cn.dat)
}
