## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

GenomeHscrSegPlot <- function(seg.dat, y.lab, y.min, y.max) {
  chr.lens <- GetChrLens(x=FALSE)
  chr.lens <- as.numeric(chr.lens)
  chr.w <- chr.lens / sum(chr.lens)
  
  plot(0, type = "n", bty = "n", xlim = c(0, 1), ylim = c(y.min, y.max),
       xlab = "Chromosome", ylab = y.lab, main = "", xaxt = "n")
  
  ww <- as.vector(rbind(chr.w, chr.w)) / 2
  chr.mids <- cumsum(ww)[(c(1:length(ww)) - 1) %% 2 == 0]
    
  lab.vals <- (c(1:length(chr.w)))
  odd.ix <- lab.vals %% 2 == 1
    
  mtext(text = lab.vals[odd.ix], side = 1, line = -0.45, at = chr.mids[odd.ix], 
        las = 1, cex = par("cex.axis") * par("cex") * 0.9)
  mtext(text = lab.vals[!odd.ix], side = 1, line = 0, at = chr.mids[!odd.ix],
        las = 1, cex = par("cex.axis") * par("cex") * 0.9)
    
  chr.offsets <- c(0, cumsum(chr.w[c(1:(length(chr.w) - 1))]), 1)
  cent.pos <- (GetCentromerePos())/sum(chr.lens) +
    chr.offsets[c(1:(length(chr.offsets) - 1))]
    
  for (i in 1:(length(chr.offsets) - 1)) {
    use.col <- ifelse(i%%2 == 1, "grey90", "white")
    rect(xleft = chr.offsets[i], ybottom = y.min, xright = chr.offsets[i + 1], 
         ytop = y.max, col = use.col, border = NA)
    lines(y = c(y.min, y.max), x = rep(cent.pos[i], 2), lty = 3, lwd = 0.5)
  }
  
  seg.cols <- GetAsSegCols(seg.dat)
  hom.color <- "darkgrey"
  het.1.color <- "red"
  het.2.color <- "blue"
  mid.color <- "darkviolet"
  pal <- colorRampPalette(c(hom.color, het.1.color, mid.color,
                            het.2.color, hom.color))
  cols <- pal(1000)
  w.d <- 0.025
    
  
  for (s in seq_len(nrow(seg.dat[["allele.segs"]]))) {
    seg.crds <- as.numeric(c(seg.dat[["allele.segs"]][s, "Start.bp"],
                             seg.dat[["allele.segs"]][s, "End.bp"]))
    chr <- as.integer(seg.dat[["allele.segs"]][s, "Chromosome"])
    genome.crds <- chr.offsets[chr] + seg.crds / chr.lens[chr] * chr.w[chr]
    means <- c(seg.dat[["allele.segs"]][s, "A1.Seg.CN"],
               seg.dat[["allele.segs"]][s, "A2.Seg.CN"], NA)
    
    het.seg.col.1 <- cols[seg.cols[s, 1]]
    het.seg.col.2 <- cols[seg.cols[s, 2]]
    
    rect(xleft = genome.crds[1], ybottom = means[1] - w.d/2, xright = genome.crds[2], 
         ytop = means[1] + w.d/2, border = NA, col = het.seg.col.1)
    rect(xleft = genome.crds[1], ybottom = means[2] - w.d/2, xright = genome.crds[2], 
         ytop = means[2] + w.d/2, border = NA, col = het.seg.col.2)
  }
}

PlotCombAnnotSegHist <- function(comb, sideways, x.max, max.q, comb.color) {
  extra <- (comb[2] - comb[1]) / 2
  
  ## FIXME: this seems like it could be simplified
  for (q in seq_along(comb)) {
    if (comb[q] > x.max) {
      break
    }
  }
  
  max.q <- q - 1
  
  for (q in seq_along(comb)) {
    ## FIXME: See above FIXME
    if (q > max.q) {
      next
    }
    
    ## ABS CN at each level
    if (!sideways) {
      abline(v=comb[q], lwd=1.5, lty=3, col=comb.color)
    } else {
      abline(h=comb[q], lwd=1.5, lty=3, col=comb.color)
    }
    
    side <- ifelse(!sideways, 3, 4)
    ## FIXME: This looks like a || instead of |
    if ((q - 1 < 10) | ((q-1) %% 2 ==1)) {
      mtext(text=(q - 1), side=side, at=comb[q], col=comb.color, line=0.2,
            cex=par("cex") * par("cex.axis"))
    }
  }
}


PlotHscrAndSeghist <- function(seg.dat, called.mode.ix, plot.hist = TRUE,
                               plot.abs.fit=FALSE, pp.calls=NA) {
  d.mar <- par("mar")
  ## eliminate right margins, small top margins
  par(mar = c(d.mar[c(1:2)], d.mar[2] * 0.25, 0))
  
  y.lab <- "Homologous copy ratio"
  y.min <- -0.05
  y.max <- 2
  
  GenomeHscrSegPlot(seg.dat, y.lab, y.min, y.max)
  abline(h = 0, lty = 3, lwd = 2)
  
  ##  seg hist plots
  x.min <- -0.05
  x.max <- 2
  bin.w <- 0.04  ## for histograms
  
  ## no left, small top margin, right 
  par(mar = c(d.mar[1], 0, d.mar[3] * 0.25, d.mar[4] * 0.5))
  
  if (plot.hist) {
    ## AS CN hist
    hom.color <- "darkgrey"
    het.1.color <- "red"
    het.2.color <- "blue"
    mid.color <- "darkviolet"
    pal <- colorRampPalette(c(hom.color, het.1.color, mid.color, het.2.color, 
                              hom.color))
    cols <- pal(1000)
    as.seg.cols <- GetAsSegCols(seg.dat)
    
    cn <- c(seg.dat[["allele.segs"]][, 6], seg.dat[["allele.segs"]][, 7])
    
    cn.ix <- (cn < x.max) & (cn > 0)
    w <- rep(seg.dat[["allele.segs"]][, 5], 2)
    w <- w / sum(w)
    use.cols <- c(as.seg.cols[, 1], as.seg.cols[, 2])
    PlotSeglenHist(cn[cn.ix], w[cn.ix], color.by = use.cols[cn.ix],
                     color.range = c(1, length(cols)), order.by = -w[cn.ix],
                     x.max = x.max, bin.w = bin.w, use.pal = cols, 
                     data.whiskers = FALSE, xlab = "", ylab = "", x.axis = FALSE,
                     xlim = c(x.min, x.max), sideways = TRUE)
    mtext(text = "Genomic fraction", line = par("mgp")[1], side = 1,
          cex = par("cex.axis"))
    abline(h = 0, lty = 3, lwd = 2)
    mtext("Summary histogram", side = 3, line = 0, adj = 0, cex = par("cex.axis"))

    if (plot.abs.fit) {
      stop("Not implemented - code as constructed when porting couldn't reach this block")
    }

  } else {
    frame()
  }
}

GetAsSegCols <- function(seg.dat) {
  if (is.null(seg.dat[["seg.expected.phase"]])) {
    seg.dat[["seg.expected.phase"]] <- matrix(0.5,
                                              nrow=nrow(seg.dat[["allele.segs"]]),
                                              ncol=2)
  }
  
  seg.col <- matrix(NA, nrow=nrow(seg.dat[["allele.segs"]]), ncol = 2)
  
  for (s in seq_len(nrow(seg.dat[["allele.segs"]]))) {
    e.state <- round((seg.dat[["seg.expected.phase"]][s, ] / 3) * 999, 0) + 1
    seg.col[s, 1] <- e.state[1]
    seg.col[s, 2] <- e.state[2]
  }
  
  return(seg.col)
}
