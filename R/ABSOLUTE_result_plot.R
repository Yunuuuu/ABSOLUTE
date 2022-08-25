## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

AbsoluteResultPlot <- function(sample.pdf.fn, seg.dat, called.mode.ix=NA,
                               verbose=FALSE) {
  pdf(sample.pdf.fn, 14, 15)
    
  ## genome HSCR segplot
  d.mar <- par("mar")
  layout(mat=matrix(1:4, nrow=2, ncol=2, byrow=TRUE),
         widths=c(6, 2), heights=c(6, 6))
  par(las=1)
  if (!is.null(seg.dat[["allele.segs"]])) {
    PlotHscrAndSeghist(seg.dat, called.mode.ix)
    frame()
    frame()
  }
  par(mar=d.mar)
    
  ## PP mode plots
  par(mfrow=c(4, 4))
  PlotModes(list(seg.dat), sample.pdf.fn, n.print=19, 
            debug.info=TRUE, add=TRUE, fig.mode=FALSE,
            sideways=FALSE, called.mode.ix=called.mode.ix, verbose=verbose)
  dev.off()
}
