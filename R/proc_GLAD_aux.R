## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

FilterSegs <- function(seg.info, min_probes=10, max_sd=100, verbose=FALSE) {  
  cn.mean <- mean(seg.info[, "copy_num"])
  cn.sd <- sd(seg.info[, "copy_num"])
  data.trunc <- cn.mean + as.numeric(max_sd) * cn.sd
  
  if (verbose) {
    print(paste(nrow(seg.info), " total segments", sep = ""))
  }
  
  rows <- (seg.info[, "n_probes"] >= min_probes) &
          (seg.info[, "copy_num"] <= data.trunc)
  
  seg.info <- seg.info[rows, ]
  
  if (verbose) {
    print(paste(nrow(seg.info), " segments remaining", sep = ""))
  }
  
  return(list(seg.info = seg.info))
}

filter_nonfinite_segments = function(seg_info, verbose=FALSE) {
  new_seg_info = seg_info[is.finite(seg_info[, "seg_sigma"]), ]
  nf_rows = nrow(seg_info) - nrow(new_seg_info)
  if ((verbose) && (nf_rows > 0)) {
    print(paste("Removed", nf_rows, "non-finite segments"))
  }
  
  return(new_seg_info)
}