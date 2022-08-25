## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

override_absolute_calls = function(segobj.list, call_override, verbose=FALSE) {
  mode_str = call_override
  mode.ix = rep(NA, length(mode_str))
  status.vec = rep(NA, length(mode_str))

  empty_ix = (call_override == "") | (is.na(call_override))
  call_override[empty_ix] = "1"
  called_ix = !is.na(as.integer(call_override))
  status.vec[called_ix ] = "called"   # PP_calls[ix,"call status"]

  mode.ix[called_ix] = as.integer(call_override[called_ix])
  mode.ix[!called_ix] = 1
  status.vec[!called_ix] = call_override[!called_ix] 

  segobj.list = SelectMatchedModes(segobj.list, mode.ix, status.vec, verbose=verbose)
  
  return( segobj.list )
}


ExtractReviewedResults = function(reviewed.pp.calls.fn, analyst.id, modes.fn, 
                                  out.dir.base, obj.name, copy_num_type, verbose=FALSE) {
  if (copy_num_type == "total") {
    set_total_funcs()
  } else if (copy_num_type == "allelic") {
    set_allelic_funcs()
  } else {
    stop("Unsupported copy number type: ", copy_num_type)
  }
  
  ## provides segobj.list
  load(modes.fn)  

  dat = read.delim(reviewed.pp.calls.fn, row.names=NULL, stringsAsFactors=FALSE, header=1,
                    check.names=FALSE)

  ## detect whether a manual override column was prepended to the front of this
  if (colnames(dat)[3] == "sample") {
     call_override = dat[, 1]
  
     pp.calls = dat[, c(3:ncol(dat))]
     rownames(pp.calls) = dat[, "sample"]
     names(call_override) = dat[, "sample"]

     found = intersect(names(segobj.list), rownames(pp.calls))
     segobj.list = segobj.list[found]
     call_override = call_override[found]

     called.segobj.list = override_absolute_calls(segobj.list, call_override)
  } else {
     if (colnames(dat)[2] != "sample") { 
       stop("Invalid reviewed.pp.calls.fn!") 
     }
     pp.calls = dat[, c(2:ncol(dat))]
     rownames(pp.calls) = dat[, "sample"]

     found = intersect(names(segobj.list), rownames(pp.calls))
     segobj.list = segobj.list[found]

     mode.ix = MatchPpModes(pp.calls, segobj.list)[["matched"]]
     called.segobj.list = ReduceSegobjListToCalled(mode.ix, segobj.list)
  }
  
  process_extract_reviewed_results(out.dir.base, called.segobj.list, pp.calls, obj.name, analyst.id)
} 

process_extract_reviewed_results = function(out.dir.base, called.segobj.list, pp.calls, 
                                            obj.name, analyst.id) {
  ## SEG_MAFs
  seg.maf.dir = file.path(out.dir.base, "reviewed", "SEG_MAF")
  dir.create(seg.maf.dir, recursive=TRUE)
  write_called_seg_maf(called.segobj.list, pp.calls, seg.maf.dir)
  
  ## PP tab
  out.fn = file.path(out.dir.base, "reviewed", paste(obj.name, ".", analyst.id, 
                                                     ".ABSOLUTE.table.txt", sep=""))
  PrintPpCallTable(called.segobj.list, out.fn)
  
  ## Called summary plot
  pdf.fn = file.path(out.dir.base, "reviewed", 
                     paste(obj.name, ".called.ABSOLUTE.plots.pdf", sep=""))
  PlotModes(called.segobj.list, pdf.fn, n.print=1)
  
  ## Called indv. RData files
   indv.called.dir = file.path(out.dir.base, "reviewed", "samples")
   dir.create(indv.called.dir, recursive=TRUE)
  
   file.base = file.path(paste(names(called.segobj.list), ".ABSOLUTE.", analyst.id, 
                               ".called", sep = ""))
   called.files= file.path(indv.called.dir, paste(file.base, "RData", sep = "."))
   for (i in seq_along(called.files)) {
      seg.obj = called.segobj.list[[i]]
      save(seg.obj, file=called.files[i])
   }
}


