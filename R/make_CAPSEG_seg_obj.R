## functions for running ABSOLUTE with total CR

total_make_seg_obj = function(dat_fn, filter_segs=FALSE, min_probes=NA, max_sd=NA, 
                              verbose=FALSE) {   
  segs_tab = read.delim(dat_fn, row.names=NULL, check.names=FALSE, stringsAsFactors=FALSE)
  seg_dat = list()
  
  nix = segs_tab[,"Chromosome"] %in% c("X", "Y", "chrX", "chrY")
  segs_tab = segs_tab[!nix, ]
  
  segs_tab[,"Chromosome"] = as.integer(gsub( "chr", "", segs_tab[,"Chromosome"]))
  segtab = segs_tab[, c("Chromosome", "Start", "End", "Num_Probes")]
  colnames(segtab) = c("Chromosome", "Start.bp", "End.bp", "n_probes")
  
  length = segtab[, "End.bp"] - segtab[, "Start.bp"]
  ## Convert from base 2 log
  copy_num = 2^(segs_tab[, "Segment_Mean"] )   
  
  ix = copy_num > 5.0
  if (verbose) {
    print(paste("Capping ", sum(ix), " segs at tCR = 5.0", sep=""))
  }
  copy_num[ix] = 5.0
  
  seg_sigma =  0.125 / sqrt(as.numeric(segs_tab[,"Num_Probes"]))   
  
  segtab = cbind(segtab, length, copy_num, seg_sigma)
  
  seg_dat$segtab = filter_nonfinite_segments(segtab, verbose=verbose)

  if (filter_segs) {
    seg_dat$segtab = FilterSegs(segtab, min_probes=min_probes, max_sd=max_sd)$seg.info
  }
  
  W = as.numeric(seg_dat$segtab[,"length"])
  W = W / sum(W)
  seg_dat$segtab = cbind(seg_dat$segtab, W)
  colnames(seg_dat$segtab)[ncol(seg_dat$segtab)] = "W"

  seg_dat$error_model = list()
  return(seg_dat)  
}

total_extract_sample_obs = function(seg.obj) {
  seg.dat = seg.obj$segtab
  d = seg.dat[, "copy_num"]
  stderr = seg.dat[, "seg_sigma"]
  W = seg.dat[,"W"]
  
  d_tx = TxData(seg.obj$error.model, d)
  seg_ix = seq_along(d)
  ## FIXME: "error.model" was originally named "HSCN_params" - double check this
  obs = list(d.tx=d_tx, d.stderr=stderr, W=W, seg.ix=seg_ix, error.model=seg.obj$error_mode, 
             segtab=seg.obj$segtab, platform="WES", data.type="TOTAL")
  
  return(obs)
}

