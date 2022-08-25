## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

CreateReviewObject = function(obj.name, absolute.files, indv.results.dir, copy_num_type, plot.modes=TRUE, verbose=FALSE)
{
  
  if (copy_num_type == "total") {
    set_total_funcs()
  } else if (copy_num_type == "allelic") {
    set_allelic_funcs()
  } else {
    stop("Unsupported copy number type: ", copy_num_type)
  }
  
  dir.create(indv.results.dir, recursive=TRUE)
  nm = file.path(indv.results.dir, paste(obj.name, ".PP-modes", sep = ""))
  modesegs.fn = paste(nm, ".data.RData", sep = "")
  pdf.fn = paste(nm, ".plots.pdf", sep = "")
  failed.pdf.fn = paste(nm, "FAILED_plots.pdf", sep = ".")
  failed.tab.fn = paste(nm, "FAILED_tab.txt", sep = ".")
  call.tab.fn = file.path(indv.results.dir,
                           paste(obj.name, "PP-calls_tab.txt", sep = "."))
  

  ## read in processed SEG / MODES results and assemble
  segobj.list = vector(mode = "list", length = length(absolute.files))
  failed.list = vector(mode="list", length=length(absolute.files))
  so.ix = 1
  fa.ix = 1
  sids = character()
  
  for (i in seq_along(absolute.files)) {
    ## read in absolute rda
    seg.out.fn = absolute.files[i]
    if (!file.exists(seg.out.fn)) {
      if (verbose) {
        cat("\n")
        print(paste("sample #", i, " result not found", sep = ""))
      }
      next
    }

    ## provides seg.dat
    load(absolute.files[i])
    SID = seg.dat$sample.name
    ## Make sure we're only using unique sample names
    if (SID %in% sids) {
      stop("Sample names must be unique")
    }
    sids = c(sids, SID)

    if (is.null(seg.dat$array.name)) {
      seg.dat$array.name = SID
    }
       
    if (is.na(seg.dat[["mode.res"]][["mode.flag"]])) {
      segobj.list[[so.ix]] = seg.dat
      names(segobj.list)[so.ix] = SID
      so.ix = so.ix + 1
      
      if (verbose) {
        cat(".")
      }
    } else {
      if (verbose) {
        print(paste("Failed list got updated", i))
      }
      failed.list[[fa.ix]] = seg.dat
      names(failed.list)[fa.ix] = SID
      failed.list[[fa.ix]][["sample.name"]] = seg.dat[["sample.name"]]
      fa.ix = fa.ix + 1
      if (verbose) {
        cat("-")
      }
    }
  }
  
  segobj.list = segobj.list[c(1:(so.ix - 1))]
  failed.list = failed.list[c(1:(fa.ix - 1))]
  
  ## sort samples by the entropy of the best solution
  mode.ent = rep(NA, length(segobj.list))
  names(mode.ent) = names(segobj.list)
  for (i in seq_along(segobj.list)) {
    mtab = segobj.list[[i]][["mode.res"]][["mode.tab"]]
    ix = which.max(mtab[, "post_LL"])
    mode.ent[i] = mtab[ix, "entropy"]
  }
  samples = names(sort(mode.ent))
  segobj.list = segobj.list[samples]
  save(segobj.list, file = modesegs.fn)


  PrintPpCallTable(segobj.list, call.tab.fn)
  
  if (plot.modes) {
    PlotModes(segobj.list, pdf.fn, n.print=3)
  }

  if (!is.null(failed.list[[1]])) {
    PlotFailedSamples(failed.list, failed.pdf.fn)
    PrintFailedTable(failed.list, failed.tab.fn)
  }
  

  return(TRUE)
}
