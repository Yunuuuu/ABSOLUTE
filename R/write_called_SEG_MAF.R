## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

write_called_seg_maf = function(called_segobj_list, pp_calls, out_dir, verbose=FALSE) {
  dir.create(out_dir) 
  
#  mode_ix = MatchPpModes(pp_calls, segobj_list)$matched
#  nix = which(is.na(mode_ix))
#  print(paste(length(nix), " samples failed PP match:", sep=""))
#  print(names(segobj_list)[nix])
#  called_segobj_list = ReduceSegobjListToCalled(mode_ix, segobj_list) 
  
  for (i in seq_along(called_segobj_list)) {
    called_segobj = called_segobj_list[[i]] 
    s_name = called_segobj$array.name
    
    at = called_segobj$mode.res$mode.tab[1, "AT"]
    abs_seg = GetAbsSegDat(called_segobj, at)
    WriteAbsSegtab(list(abs_seg), s_name, 
                   file.path(out_dir, paste(s_name, "segtab.txt", sep=".")))
    ## MAF
    if (!is.null(called_segobj$mode.res$muts.post.prs)) {
      maf_out_fn = file.path(out_dir, paste(s_name, "_ABS_MAF.txt", sep=""))
            
      modeled = called_segobj$mode.res$muts.post.prs[, , 1, drop=TRUE] 
      ## check case for only 1 mut
      if (is.null(dim(modeled))) { 
        modeled = matrix(modeled, nrow=1)
        colnames(modeled) = dimnames(called_segobj$mode.res$muts.post.prs)[[2]]  
      }
      mut_dat = cbind(called_segobj$mut.cn.dat, modeled)      
      
      ## SNVs rescaling to cell-fraction
      ccf_grid = seq(0, 1, by=0.01)
      alpha = mut_dat[1, "purity"]
      ccf_dens = matrix(NA, nrow=nrow(mut_dat), ncol=length(ccf_grid))
      for (i in seq_along(nrow(mut_dat))) {
        if (mut_dat[i, "q_hat"] == 0) { 
          next 
        }
        ccf_dens[i, ] = calc_ccf_posterior_grid(mut_dat[i, "alt"], mut_dat[i, "ref"], 
                                                alpha, mut_dat[i, "q_hat"], ccf_grid)
      }
      colnames(ccf_dens) = paste("CCF_", ccf_grid, sep="")
      out_mut_dat = cbind(mut_dat, ccf_dens)
            
      write.table(file=maf_out_fn, out_mut_dat, row.names=FALSE, sep="\t", quote=FALSE)
    }
    if (verbose) {   
      cat(".")
    }
  }
}

WriteAbsSegtab <- function(seg, s_name, out_fn) {
  for (s in seq_along(seg)) {
    sample <- rep(s_name, nrow(seg[[s]]))
    s_tab <- cbind(sample, seg[[s]])
    
    ## colames only for 1st sample
    app <- s > 1   
    col <- s == 1
    
    write.table(s_tab, file=out_fn, col.names=col, append=app, row.names=FALSE,
                quote=FALSE, sep="\t")
  }   
}

