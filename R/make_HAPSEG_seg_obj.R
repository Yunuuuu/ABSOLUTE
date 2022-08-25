## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

AllelicMakeSegObj <- function(seg.dat.fn, filter_segs=FALSE, min_probes=NA, max_sd=NA,
                              use.gcm=FALSE, verbose=FALSE) {   
  ## provides seg.dat
  if (!file.exists(seg.dat.fn)) {
    stop("seg.dat.fn does not exist")
  }
  load(seg.dat.fn)
  if (!exists("seg.dat")) {
    stop("Invalid object contained in seg.dat.fn")
  }

  ## Detect pre-1.0 HAPSEG inputs and convert
  if ("EM_res" %in% names(seg.dat)) {
    if (verbose) {
      print("Detected a pre-1.0 HAPSEG input, converting to 1.0 format ....")
    }
    seg.dat = RenameSegDatFields(seg.dat)
  }
  
  ## Replace GCM CN segs with sum of AS segs
  if (use.gcm == FALSE) {
    new.seg.info <- seg.dat[["allele.segs"]][, c("Chromosome",
                                                 "Start.bp",
                                                 "End.bp", "n_probes",
                                                 "length", "A1.Seg.CN",
                                                 "tCN.Seg.sd")]
    new.seg.info[, "A1.Seg.CN"] <- apply(seg.dat[["allele.segs"]][, c("A1.Seg.CN", "A2.Seg.CN")], 1, sum) / 2
    colnames(new.seg.info)[c(6, 7)] <- c("copy_num", "seg_sigma")
    seg.info <- new.seg.info
  } else {
    seg.info <- seg.dat[["seg.info"]]
  }
  seg.dat[["seg.info"]] <- filter_nonfinite_segments(seg.info, verbose=verbose)
  if (filter_segs == TRUE) {
    seg.dat[["seg.info"]] <- FilterSegs(seg.dat[["seg.info"]], min_probes, max_sd,
                                        verbose=verbose)[["seg.info"]]
  }
  
  W <- as.numeric(seg.dat[["seg.info"]][, "length"])
  W <- W / sum(W)
  seg.dat[["seg.info"]] <- cbind(seg.dat[["seg.info"]], W)
  colnames(seg.dat[["seg.info"]])[ncol(seg.dat[["seg.info"]])] <- "W"
  
  ## also create AS seg.info
  as.1 <- seg.dat[["allele.segs"]][, c("Chromosome", "Start.bp", "End.bp",
                                       "n_probes", "length", "A1.Seg.CN",
                                       "AS.Seg.sd")]
  as.2 <- seg.dat[["allele.segs"]][, c("Chromosome", "Start.bp", "End.bp",
                                       "n_probes", "length", "A2.Seg.CN",
                                       "AS.Seg.sd")]
  
  colnames(as.1)[6] <- "copy_num"
  colnames(as.2)[6] <- "copy_num"
  
  ## remember each segment's mate
  seg.ix <- c(1:nrow(as.1))
  as.1 <- cbind(as.1, seg.ix)
  as.2 <- cbind(as.2, seg.ix)
  
  as.seg.dat <- rbind(as.1, as.2)
  colnames(as.seg.dat)[7] <- "seg_sigma"
  colnames(as.seg.dat)[8] <- "seg.ix"
  seg.dat[["as.seg.dat"]] <- filter_nonfinite_segments(as.seg.dat, verbose=verbose)
  
  if (filter_segs) {
    seg.dat$as.seg.dat = FilterSegs(seg.dat$as.seg.dat, min_probes=min_probes, max_sd=max_sd,
                                verbose=verbose)$seg.info
  }
  
  W <- as.numeric(seg.dat[["as.seg.dat"]][, "length"])
  W <- W / sum(W)
  seg.dat[["as.seg.dat"]] <- cbind(seg.dat[["as.seg.dat"]], W)
  colnames(seg.dat[["as.seg.dat"]])[ncol(seg.dat[["as.seg.dat"]])] <- "W"
  
  ## add error-model parameters from EM-fit
  em <- seg.dat[["em.res"]]
  
  if (is.null(seg.dat[["em.res"]][["theta"]])) {
    ## This indicates an older hapseg output, these values aren't
    ## contained by theta
    seg.dat[["error.model"]] <- list()
    seg.dat[["error.model"]][["sigma.eta"]] <- em[["sigma.eta"]]
    seg.dat[["error.model"]][["sigma.nu"]] <- em[["sigma.nu"]]
    seg.dat[["error.model"]][["het.cov"]] <- em[["het.cov"]]
    seg.dat[["error.model"]][["nu"]] <- em[["nu"]]
    seg.dat[["error.model"]][["AT"]] <- em[["at"]]
    seg.dat[["error.model"]][["BG"]]<- em[["bg"]]
  } else {
    seg.dat[["error.model"]] <- seg.dat[["em.res"]][["theta"]]
  }
  
  seg.dat[["error.model"]][["loglik"]] <- em[["theta"]][["loglik"]]

  ## Remove multiple fields from the resulting list as we no longer
  ## need them
  remove.fields <- c("em.res", "c.hat.1", "g.0",
                     "c.hat.phased.unfixed", "h.switch.ix",
                     "merged.loci", "final.merge.prob")
  seg.dat <- seg.dat[setdiff(names(seg.dat), remove.fields)]
  
  return(seg.dat)
}

AllelicExtractSampleObs <- function(seg.obj) {
  seg.dat <- seg.obj[["as.seg.dat"]]
  d <- seg.dat[, "copy_num"]
  stderr <- seg.dat[, "seg_sigma"]
  W <- seg.dat[, "W"]
  if("seg.ix" %in% colnames(seg.dat)) {
    seg.ix <- seg.dat[,"seg.ix"]
  } else {
    seg.ix <- NA
  }
  
  hscn.params <- seg.obj[["error.model"]]
  
  ## initialize AT result to be fit
  ## FIXME: This appears to be another spot w/ name issues between
  ## versions
  if ("AT" %in% names(hscn.params)) {
    hscn.params[["at"]] <- hscn.params[["AT"]]
  }
  hscn.params[["fit.at"]] <- hscn.params[["at"]] 
  
  d.tx <- TxData(hscn.params, d)
  
  obs <- list(d.tx=d.tx, d.stderr=stderr, W=W, seg.ix=seg.ix,
              error.model=hscn.params, segtab=seg.obj[["as.seg.dat"]],
              platform="ARRAY", data.type="ALLELIC")
  
  return(obs)
}





