## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

MatchPpModes = function(pp.calls, segobj_list) {
  N = length(segobj_list)
  mode.ix = rep(NA, N)
  wrong.ix = c()
  names(mode.ix) = names(segobj_list)
  
  for (i in seq_len(N)) {
    res = MatchPpMode(pp.calls, names(segobj_list)[i],
                       segobj_list[[i]])
    mode.ix[i] = res[1]

    if (!is.na(res) && res[2]) {
      wrong.ix = c(wrong.ix, i)
    }
  }
  
  return(list(matched=mode.ix, wrong.ix=wrong.ix))
}

MatchPpMode = function(pp.calls, sid, seg.dat) {
  if ((any(is.na(pp.calls[sid, c("purity", "ploidy")]))) ||
      (!is.na(seg.dat[["mode.res"]][["mode.flag"]]))) {
    return(NA)
  }
  
  mode.tab = seg.dat[["mode.res"]][["mode.tab"]]
  vals = mode.tab[, c("alpha", "genome mass"), drop=FALSE]   
  
  call.vals = pp.calls[sid, c("purity", "ploidy")]
  call.vals = matrix(as.numeric(call.vals), ncol=2, nrow=nrow(mode.tab), byrow=TRUE)
  
  min.mode = which.min(rowSums((vals - call.vals)^2))
  
  if ((abs(pp.calls[sid, "ploidy"] - vals[min.mode, "genome mass"]) > 1) ||
      (abs(pp.calls[sid, "purity"] - vals[min.mode, "alpha"]) > 0.1)) {
    return(NA)  
  }

  mode.ix = min.mode       
  
  ## FIXME: Scott requested the old variant stays here. Note that 'wrong'
  ## isn't really used anywhere
  ## if (mode.ix != which.max(mode.tab[, "dens"])) {
  ##    wrong = TRUE
  ##  } else {
  ##    wrong = FALSE
  ##  }
 wrong = TRUE
  
  return(c(mode.ix, wrong))
}

ReduceSegobjListToCalled = function(mode.ix, segobj_list) {
  segobj_list = segobj_list[!is.na(mode.ix)]
  mode.ix = mode.ix[!is.na(mode.ix)]
  
  for (i in seq_along(segobj_list))  {
    segobj_list[[i]][["mode.res"]] = ReorderModeRes(segobj_list[[i]][["mode.res"]],
                                                     mode.ix[i])
  }
  
  return(segobj_list)
}

SelectMatchedModes = function(segobj_list, mode.ix, status.vec, verbose=FALSE) {
   new.segobj.list = segobj_list[!is.na(mode.ix)]
   mode.ix = mode.ix[!is.na(mode.ix)]

   for (i in seq_along(mode.ix)) {
      segobj = new.segobj.list[[i]]
      segobj[["mode.res"]][["call.status"]] = status.vec[i]

      if (!(status.vec[i] %in% c("non-clonal", "non-aneuploid",
                                "low purity", "FAILED"))) {
         mode.res = segobj[["mode.res"]]
         ix = mode.ix[i]
         segobj[["mode.res"]] = ReorderModeRes(mode.res, ix)
      }

      new.segobj.list[[i]] = segobj
      if (verbose) {
        cat(".")
      }
   }
   if (verbose) {
     cat("\n")
   }
   
   return(new.segobj.list)
}
