## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

RunAbsolute = function(seg.dat.fn, sigma.p, max.sigma.h, min.ploidy, max.ploidy,
                        primary.disease, platform, sample.name, results.dir,
                        max.as.seg.count, max.non.clonal, max.neg.genome,
                        copy_num_type, maf.fn = NULL, min.mut.af = NULL,
                        output.fn.base=NULL, verbose = FALSE) {    
  ## Other constants
  kQ = 15
  kTauDom = c(min.ploidy, max.ploidy)
  kSigmaHDom = c(0, max.sigma.h)
  kPiSomThetaQ = c(100, 50, rep(2, (kQ - 2)))
  mut_class_w = list(SM = 0.5, GL = 0, SC = 0.5, OL = 1e-3, Pi_SM = 15, Pi_SC = 15)

  platform = match.arg(platform, c("SNP_6.0", "Illumina_WES", "SNP_250K_STY"))
  if (platform %in% c("SNP_6.0", "SNP_250K_STY")) {
    filter_segs = FALSE
  } else if (platform == "Illumina_WES") {
    filter_segs = TRUE
  } else {
    stop("Unsupported platform: ", platform)
  }
  min_probes = 10
  max_sd = 100
  
  copy_num_type = match.arg(copy_num_type, c("allelic", "total"))
  if (copy_num_type == "total") {
    pi_theta_qz = list(W = c(1, 25, 100, 25, 10, 5, rep(1, kQ - 6), 10), PC = 0.05)
    set_total_funcs()
  } else if (copy_num_type == "allelic") {
    pi_theta_qz = list(W = c(25, 100, 25, 10, 5, rep(1, kQ - 5), 10), PC = 0.05)
    set_allelic_funcs()
  } else {
    stop("Unsupported copy number type: ", copy_num_type)
  }
  
  ## FIXME: I think this should be exposed?
  SubclonalPost <<- ABSOLUTE:::UnifSubclonalPost
  
  ## sigma.p becomes sigma.h
  ## FIXME: Should this be changed throughout?
  sigma.h = sigma.p
  
  tmp.dir = file.path(results.dir, "tmp")

  dir.create(tmp.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(results.dir, recursive=TRUE, showWarnings=FALSE)
  
  seg.dat = MakeSegObj(seg.dat.fn, min_probes=min_probes, max_sd=max_sd, 
                        filter_segs=filter_segs, verbose=verbose)
  seg.dat[["primary.disease"]] = primary.disease
  seg.dat[["group"]] = DetermineGroup(primary.disease)
  seg.dat[["platform"]] = platform
  seg.dat[["sample.name"]] = as.character(sample.name)
  if (is.null(seg.dat$array.name)) {
    seg.dat$array.name = seg.dat$sample.name
  }
  seg.dat[["maf.fn"]] = maf.fn
  
  ## either allelic or total CR, to be modeled.
  seg.dat[["obs.scna"]] = ExtractSampleObs(seg.dat)
  
  ## check for QC failure modes
  seg = seg.dat[["obs.scna"]][["segtab"]]
  ## expected copy-number, should be 1.0
  e.cr = sum(seg[, "W"] * seg[, "copy_num"])
  if (verbose) {
    print(paste("Expected copy-ratio = ", round(e.cr, 5), sep=""))
  }
  
  mode.res = list(mode.flag = NA)
  
  if (nrow(seg) > max.as.seg.count) {
    mode.res[["mode.flag"]] = "OVERSEG"
  }
  
  if ((e.cr < 0.75) || (e.cr > 1.25)) {
    mode.res[["mode.flag"]] = "E_CR_SCALE"
  }
  
  if (is.na(mode.res[["mode.flag"]])) {
    ## check for MAF describing somatic mutations
    maf = NULL
    if ((!is.null(maf.fn)) && (file.exists(maf.fn))) {
      maf = read.delim(maf.fn, row.names = NULL, stringsAsFactors = FALSE, 
                        check.names = FALSE, na.strings = c("NA", "---"),
                        blank.lines.skip=TRUE, comment.char="#")
    } else {
      if (verbose) {
        print(paste("MAF file: ", maf.fn, " not found.", sep = ""))
      }
    }
    
    ## find initial purity/ploidy solutions for debugger
    data(ChrArmsDat, package = "ABSOLUTE")
    if ((!is.null(maf)) && (nrow(maf) > 0)) {
      if (is.null(min.mut.af)) {
        stop("min.mut.af must be defined")
      }
      mut.cn.dat = CreateMutCnDat(maf, seg.dat, min.mut.af, verbose=verbose)
      mut = list(mut_CN_dat = mut.cn.dat, Pi_som_theta_Q = kPiSomThetaQ,
                 mut_class_W = mut_class_w)      
    } else {
      mut = NA
    }

    mode.res = FitSample(seg.dat, mut, kQ, pi_theta_qz, sigma.h, kTauDom,
                         kSigmaHDom, chr.arms.dat, verbose=verbose)
    mode.res = apply_subclonal_scna_model(seg.dat, mode.res, verbose=verbose)
    
    if (inherits(mode.res, "try-error")) {
      mode.res = list(mode.flag = "FAIL")
    }
  }

  if (is.na(mode.res[["mode.flag"]])) {
    bad.ix = GenomeHetFilter(seg.dat[["obs.scna"]], mode.res, max.non.clonal,
                              max.neg.genome, kQ, verbose=verbose)
    if (sum(bad.ix) == nrow(mode.res[["mode.tab"]])) {
      mode.res = list(mode.flag="ALPHA_TAU_DOM")
    } else {
      mode.res = ReorderModeRes(mode.res, !bad.ix)
    }
  }
  
  if (is.na(mode.res[["mode.flag"]])) {
    ## implementation below is order-dependent.
    ##  Will break if no Kar model but MAF is present
    ## 1 - apply karyotype model
    data(ChrArmPriorDb, package="ABSOLUTE")
    model.id = ifelse(seg.dat[["group"]] %in% names(train.obj),
                       seg.dat[["group"]], "Primary")
    mode.res = ApplyChrArmPrior(mode.res, model.id, train.obj,
                                 verbose=verbose)
    
    ## 2 - apply mutation model
    if ((!is.null(maf)) && (nrow(maf) > 0)) {
      if (is.null(min.mut.af)) {
        stop("maf was not NULL, but min.mut.af was")
      }
      seg.dat[["mut.cn.dat"]] = mut.cn.dat
      mode.res = ApplySomaticMutsModel(mode.res, seg.dat$obs_SCNA, mut.cn.dat,
                                        kPiSomThetaQ, mut_class_w, kQ, verbose=verbose)
    }
    
    mode.res = WeighSampleModes(mode.res)
    mode.res[["call.status"]] = GetCallStatus(mode.res,
                                               seg.dat[["obs.scna"]][["W"]])
  }
  
  seg.dat[["mode.res"]] = mode.res

  if (is.null(output.fn.base)) {
    output.fn.base = ifelse(is.null(seg.dat$array.name), sample.name, seg.dat$array.name)
  }
    
  file.base = paste(output.fn.base, ".ABSOLUTE", sep = "")
  
  if (is.na(mode.res[["mode.flag"]])) {
    sample.pdf.fn = file.path(results.dir,
                               paste(file.base, "plot.pdf", sep = "_"))
    AbsoluteResultPlot(sample.pdf.fn, seg.dat)
  } else {
    if (verbose) {
      print("Mode flag is NA, not generating plots. Sample has failed ABSOLUTE")
    }    
  }
  
  seg.dat$version = 1.1
  
  save(seg.dat, file = file.path(results.dir, paste(file.base, "RData", sep = ".")))
  
  return(TRUE)
}

RenameSegDatFields = function(seg.dat) {
  ## FIXME: Total hack, used to go back and forth between the old
  ## hapseg & current absolute
  
  names(seg.dat) = ConvertNameStyle(names(seg.dat))
  names(seg.dat[["em.res"]]) = ConvertNameStyle(names(seg.dat[["em.res"]]))
 # names(seg.dat[["obs.scna"]]) = ConvertNameStyle(names(seg.dat[["obs.scna"]]))
 # names(seg.dat[["error.model"]]) = ConvertNameStyle(names(seg.dat[["error.model"]]))
  #names(seg.dat[["mode.res"]]) = ConvertNameStyle(names(seg.dat[["mode.res"]]))
  
  return(seg.dat)
}

ConvertNameStyle = function(names) {
  return(tolower(gsub("_", ".", names)))
}
