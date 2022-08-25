## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

GetModeColors <- function() {
  return(c("springgreen3", "dodgerblue", "maroon1", "goldenrod1",
           "orangered2", 4, 6, "salmon", "yellowgreen", "cyan",
           "sienna", "rosybrown", rep(1,100) ))
}

PlotModes <- function(segobj.list, pdf.fn, 
                      n.print = NA, write.tab=TRUE, debug.info = TRUE,
                      add = FALSE, fig.mode = FALSE, sideways = FALSE, 
                      low.w=FALSE, called.mode.ix=NA, verbose=FALSE) {
  pal <- colorRampPalette(c("chocolate4", "cyan4"))
  
  if (fig.mode) {
    add <- TRUE
    debug.info <- FALSE
  }
  
  Q <- 15
  alpha.dom <- c(0, 1)
  tau.dom <- c(1, 10)
  mode.colors <- GetModeColors()
  x.max <- 2.25
  binW <- 0.025
  
  if (!add) {
    ## Multi-sample summary mode
    pdf(pdf.fn, 3.5 * (n.print + 2), 15)
    par(mfrow = c(5, n.print + 3))
  }
  
  for (s in seq_len(length(segobj.list))) {
    mode.tab <- segobj.list[[s]][["mode.res"]][["mode.tab"]]
    seg.z.tab <- segobj.list[[s]][["mode.res"]][["seg.z.tab"]]
    
    if (nrow(mode.tab) == 0) {
      for (i in seq_len(n.print + 2)) {
        frame()
      }
      next
    }
    
    obs <- ExtractSampleObs(segobj.list[[s]])
    SN <- segobj.list[[s]][["sample.name"]]
    
    if (is.na(n.print)) {
      s.n.print <- nrow(mode.tab)
    } else {
      s.n.print <- min(nrow(mode.tab), n.print)
    }
    
    ## Note: In Scott's code this was initialized to 0 and then
    ## he ran a for loop 1:s.n.print, incrementing s.mix. In case
    ## I wonder in the future why this is so different ....
    s.mix <- s.n.print
    
    if (!fig.mode) {
      model.id <- segobj.list[[s]][["group"]]
      
      PostPlot(mode.tab, mode.colors, alpha.dom, tau.dom, SN, obs,
               called.mode.ix, debug.info,
               call.status=segobj.list[[s]][["mode.res"]][["call.status"]],
               model.id=model.id)
      PpModeScorerBarplot(mode.tab, mode.colors, obs, n.print)
      
      if (length(segobj.list) == 1) {
        frame()
        frame()
      }
    }
    
    n.plot <- min(n.print, NROW(mode.tab))
    mode.tab <- mode.tab[c(1:n.plot), , drop = FALSE]
    
    for (i in seq_len(n.plot)) {
      tau <- mode.tab[i, "tau"]
      alpha <- mode.tab[i, "alpha"]
      delta <- alpha / (2 * (1 - alpha) + alpha * tau)
      b <- 2 * (1 - alpha) / (2 * (1 - alpha) + alpha * tau)
      mode.info <- mode.tab[i, ]
      
      obs[["error.model"]][["fit.at"]] <- mode.tab[i, "AT"]
      comb <- InvTxData(obs[["error.model"]],
                        GetCopyRatioComb(Q, delta, b, obs[["error.model"]]))
      seg.dat <- InvTxData(obs[["error.model"]], obs[["d.tx"]])
      
      last.plot <- ifelse(i == n.plot, TRUE, FALSE)
      mid.plot <- ifelse(i == ceiling(n.plot / 2), TRUE, FALSE)
      first.plot <- ifelse(i == 1, TRUE, FALSE)
      
      ## layout works for multi-sample mode
      if ((length(segobj.list) > 1) && (!is.null(segobj.list[[s]][["mut.cn.dat"]])) && 
          (i > 1)) {
        frame()
        frame()
      }
      
      if (!is.null(segobj.list[[s]][["mode.res"]][["ab.tab"]])) {
        comb.ab <- segobj.list[[s]][["mode.res"]][["ab.tab"]][i, ]
      } else {
        comb.ab <- NA
      }
      
      ModeCombPlot(seg.dat, obs[["W"]], obs[["d.stderr"]], mode.info,
                   seg.z.tab[i,], comb.ab, mode.colors[i], comb, last.plot,
                   mid.plot, first.plot, x.max, debug.info, fig.mode, sideways,
                   pal)
      
      if (!is.na(called.mode.ix)) {
        if (i == called.mode.ix) {
          ## highlight called solution            
          mtext(text = "*", side = 3, at = 0, col = "red", line = -1, cex = 3 * 
                par("cex") * par("cex.axis"))
        }
      }
      
      if (!is.null(segobj.list[[s]][["mut.cn.dat"]])) {
        mut.cn.dat <- segobj.list[[s]][["mut.cn.dat"]]
        modeled <- segobj.list[[s]][["mode.res"]][["muts.post.prs"]][, , i, drop=TRUE]
        ## check case for only 1 mut
        if (is.null(dim(modeled))) {
          modeled <- matrix(modeled, nrow=1)
          colnames(modeled) <- dimnames(segobj.list[[s]][["mode.res"]][["muts.post_prs"]])[[2]]
        }
        modeled.mut.dat <- cbind(mut.cn.dat, modeled)
        PlotSomaticMutDensities(modeled.mut.dat, segobj.list[[s]], 1,
                                mode.colors[i], min.cov = 3, verbose=verbose)
      }
    }
    if ((length(segobj.list) > 1) &&
        (is.null(segobj.list[[s]][["mut.cn.dat"]]))) {
      if (i < n.print) {
        for (j in seq_len(n.print - i)) {
          frame()
        }
      } 
        frame()
      
    }
  }
  if (!add) {
    dev.off()
  }
}

PlotCalledPpSolution <- function(called.segobj, log.w=FALSE) {
   mode.tab <- called.segobj[["mode.res"]][["mode.tab"]]
   
   tau <- mode.tab[1, "tau" ]
   alpha <- mode.tab[1, "alpha"]
   delta <- alpha / (2 * (1 - alpha) + alpha * tau) 
   b <- 2 * (1 - alpha) / (2 * (1 - alpha) + alpha * tau) 

   mode.info <- mode.tab[1, ]
   kXMax <- 2.25

   seg.z.tab <- called.segobj[["mode.res"]][["seg.z.tab"]][1, ]

   if (!is.null(called.segobj[["mode.res"]][["ab.tab"]])) {
     comb_AB <- called.segobj[["mode.res"]][["ab.tab"]][1, ]
   } else{
     comb_AB <- NA
   }
   
   obs <- ExtractSampleObs(called.segobj)

   obs[["error.model"]][["fit.at"]] <- mode.tab[1, "AT"]
   comb <- InvTxData(obs[["error.model"]],
                     GetCopyRatioComb(Q, delta, b, obs[["error.model"]]))
   seg.dat <- InvTxData(obs[["error.model"]], obs[["d.tx"]])

   ModeCombPlot(seg.dat, obs[["W"]], obs[["d.stderr"]],
                mode.info, seg.z.tab, comb_AB, 1, comb,
                FALSE,FALSE, FALSE, x.max=kXMax,
                debug.info=TRUE, fig.mode=FALSE, sideways=FALSE, pal=pal)
}

PostPlot <- function(mode.tab, mode.colors, alpha.dom, tau.dom, sample.name, 
                      obs, called.mode.ix, debug.info, call.status = NA,
                      model.id = NA) {
  n.print <- nrow(mode.tab)
  
  mode.tab <- mode.tab[c(1:n.print), , drop = FALSE]
  mode.colors <- mode.colors[c(1:NROW(mode.tab))]
    
  a1 <- mode.tab[1, "alpha"]
  t1 <- mode.tab[1, "tau"]
  den <- 2 * (1 - a1) + a1 * t1
  delta <- a1 / den
  b <- 2 * (1 - a1) / den
    
  delta.set <- c(delta / 2, delta, 2 * delta)
  b.set <- c(b, b - 2 * delta, b + 2 * delta)
  b.set <- b.set[b.set >= 0 & b.set < 1]
    
  old.mar <- par("mar")
  if (debug.info) {
    ## add margin line to top
    par(mar = old.mar + c(0, 0, 1, 0))  
  }
  
  ylab <- substitute( paste( "Fraction cancer nuclei", (hat(alpha)), sep="") )
  xlab <- substitute( paste( "Ploidy", (hat(tau)), sep="") )

  plot(0, type = "n", ylab = ylab, xlab = xlab,
       xlim = tau.dom, ylim = alpha.dom, main = "", bty = "n", las = 1,
       cex = 1.5, xaxt = "n")
  
  at <- seq(2, max(tau.dom), by = 2)
  axis(side = 1, at = at, labels = paste(at, "N", sep = ""))
  
  for (i in seq_along(b.set)) {
    curve(AlphaBFunc(x, b.set[i]), from = tau.dom[1], to = tau.dom[2],
          col = "grey30", add = TRUE, lty = 3)
  }
  
  for (i in seq_along(delta.set)) {
    curve(AlphaFunc(x, delta.set[i]), from = tau.dom[1], to = tau.dom[2],
          col = "grey30", add = TRUE, lty = 3, n = 500)
  }
  
  points(mode.tab[, "tau"], mode.tab[, "alpha"], col = "black", bg = mode.colors, 
         pch = 21, cex = par("cex") * 2.0)
  
  if (!is.na(called.mode.ix)) {
    ix <- called.mode.ix
    points(mode.tab[ix, "tau"], mode.tab[ix, "alpha"], col = "red", pch = "*", 
           cex = par("cex") * 4)
  }
  
  if (debug.info) {
    title(sample.name, line = 3)
    mtext(call.status, line = 0, adj = 0)
    mtext(model.id, line = 0, adj = 1)
  }
  
  if ((obs[["platform"]] == "ARRAY") && debug.info) {
    hscn.params <- obs[["error.model"]]
    mtext(substitute(paste(sigma[eta] == x, ",", sigma[epsilon] == y, sep = ""), 
                     list(x = round(hscn.params[["sigma.eta"]], 3),
                          y = round(hscn.params[["sigma.nu"]], 3))),
          line = 1.5, adj = 0)
  }
  
  par(mar = old.mar)
}


PpModeScorerBarplot <- function(mode.tab, mode.colors, obs, n.print) {
  n.print <- min(n.print, nrow(mode.tab))
  mode.colors <- mode.colors[c(1:NROW(mode.tab))]
  
  if (("somatic.mut.ll" %in% colnames(mode.tab)) &&
      (!any(is.na(mode.tab[1:n.print, "somatic.mut.ll"])))) {
    mat <- mode.tab[1:n.print, c("log.partition", "clust_LL", "somatic.mut.ll"), 
                    drop = FALSE]
    colnames(mat) <- c("SCNAs", "karyotype", "SSNVs")
  } else {
    mat <- mode.tab[c(1:n.print), c("log.partition", "clust_LL"), drop = FALSE]
    colnames(mat) <- c("SCNAs", "karyotype")
  }
  
  ix <- apply(!is.finite(mat), 1, sum) > 0
  mat <- mat[!ix, , drop=FALSE]
  
#  mat <- mat - min(mat)
  mat <- t( t(mat) - apply(mat, 2, min) )
  mat <- mat + max(mat) * 0.1
  mat <- cbind(mat, rowMeans(mat))
  colnames(mat)[ncol(mat)] <- "combined"
  
  barplot(mat, beside = TRUE, col = mode.colors[1:n.print[!ix]], axes = FALSE, 
          ylab = "", space = c(0, 2), cex.names = par("cex.axis"))
  
  mtext("Log-likelihood", side = 2, line = 1, las = 3,
        cex = par("cex") * par("cex.axis"))
  mtext("Model-based evaluation", side = 3, line = 0, adj = 0,
        cex = par("cex.axis"))
  
  axis(side=2, labels=FALSE)
}

ModeCombPlot <- function(d, W, seg.stderr, mode.info, seg.z, comb.ab, 
                         mode.color, comb, last.plot, mid.plot, first.plot,
                         x.max, debug.info, fig.mode, sideways, pal) {
  ix <- (d >= 0) & (d < x.max)
  colpal <- pal(1000)

  ## plot modes in vertical column
  col <- FALSE
  old.mar <- par("mar")

  xlab="Copy ratio"
  
  if (debug.info) {
    par(mar = old.mar + c(0, 0, 1, 0))  ## add margin line to top
  }
  
  if (fig.mode) {
    x.max <- 2
    ix <- (d >= 0) & (d < x.max)
    
    if (col) {
      x.ax.labs <- ifelse(last.plot, TRUE, FALSE)
    } else {
      x.ax.labs <- TRUE
    }
    
    if (col) {
      xlab <- ifelse(last.plot, "Allelic copy ratio", "")
    } else {
      xlab <- ifelse(first.plot, "Allelic copy ratio", "")
    }
    
    if (col) {
      ylab <- "Genomic fraction"
    } else {
      ylab <- ifelse(mid.plot, "Genomic fraction", "")
    }
    
    if (!col) {
      y.ax.labs <- first.plot
    } else {
      y.ax.labs <- last.plot
    }
    
    if (sideways) {
      tmp <- xlab
      xlab <- ylab
      ylab <- tmp
    }
    if (sideways) {
      tmp <- x.ax.labs
      x.ax.labs <- y.ax.labs
      y.ax.labs <- tmp
    }
    
    if (!fig.mode) {
      order.by <- W[ix]
    } else {
      order.by <- seg.z[ix]
    }
    
    PlotSeglenHist(d[ix], W[ix], color.by = seg.z[ix], order.by = order.by, 
                   color.range = c(0, 1), use.pal = colpal, bin.w = 0.04,
                   x.max = x.max, data.whiskers = debug.info, 
                   ylab = ylab, xlab = xlab, x.ax.labs = x.ax.labs,
                   y.ax.labs = y.ax.labs, 
                   sideways = sideways, border = !fig.mode)
    
    if (!col) {
      if (first.plot) {
        mtext("Candidate interpretations of copy profile", side = 3, line = 0, 
              adj = 0, cex = par("cex.axis"))
      }
    }
  } else {
    PlotSeglenHist(d[ix], W[ix], color.by = seg.z[ix], color.range = c(0, 1), 
                   use.pal = colpal, bin.w = 0.025, x.max = x.max,
                   data.whiskers = FALSE, xlab = xlab, 
                   sideways = sideways)
  }
  
  msg1 <- substitute( paste( purity(hat(alpha)) == x1, ", ",
			   "ploidy: ", hat(tau) == x2, ", ",
                           hat(tau)[g] == x3,  sep = ""),
                     list(
                          x1 = round(mode.info["alpha"], 2),
                          x2 = round(mode.info["tau"], 2),
                          x3 = round(mode.info["genome mass"], 2)))
  

#print(msg1)


  msg2 <- substitute(paste( hat(sigma[H]) == x1, ", ",
   			    hat(theta[Z]) == x2, ", ", "%", Het == x3, sep = ""), 
                     list( x1 = round(mode.info["sigma.h.hat"], 3),
			   x2 = round(mode.info["theta.z.hat"], 2),
                           x3 = round(100 * mode.info["frac.het"], 2)))
  
#print(msg2)

  msg3 <- substitute(paste("SCNAs" == x6, ", ", "Kar" == x7, ", ", "SSNVs" == x8, 
                           sep = ""),
                     list(x6 = round(mode.info["log.partition"], 2),
                          x7 = round(mode.info["clust_LL"], 2),
                          x8 = round(mode.info["somatic.mut.ll"], 2)))
  
#print(msg3)

  if (debug.info) {
    mtext(msg1, line = 4, cex = 0.75, adj = 0)
    mtext(msg2, line = 3, cex = 0.75, adj = 0)
    mtext(msg3, line = 2, cex = 0.75, adj = 0)
  }
  
  extra <- (comb[2] - comb[1]) / 2
  
  for (q in seq_along(comb)) {
    if (comb[q] > x.max) {
      break
    }
    
    ## calculate frac genome to the left of comb[q]
    q_ix <- ifelse(q > 1, q - 1, q)
    seg.ix <- d < comb[q_ix] - extra
    cum_gen <- sum(W[seg.ix])
    if (cum_gen > 0.99) {
      break
    }
  }
  max.q <- q - 1
  top.axis <- FALSE
  if (top.axis) {
    axis(side = 3, at = comb[c(1:max.q)], labels = c(1:max.q) - 1,
         col = mode.color, col.ticks = mode.color)
  }
  
  for (q in seq_along(comb)) {
    if (q > max.q) {
      next
    }
    
    ## ABS CN at each level
    if (!sideways) {
      abline(v = comb[q], lwd = 1.5, lty = 3, col = mode.color)
    } else {
      abline(h = comb[q], lwd = 1.5, lty = 3, col = mode.color)
    }
    
    side <- ifelse(!sideways, 3, 4)
    
    if (!top.axis) {
      if (q - 1 < 10 | (q - 1)%%2 == 1) {
        mtext(text = (q - 1), side = side, at = comb[q], col = mode.color, 
              line = 0.2, cex = par("cex") * par("cex.axis"))
      }
    }
    
    ## % AB in each level
    if (debug.info & !is.na(comb.ab)) {
      q.ab <- round(comb.ab[q], 2)
      if (is.finite(q.ab) & q.ab > 0) {
        mtext(text = paste("(", q.ab, ")", sep = ""), at = comb[q], col = 1, 
              line = 1, cex = 0.5)
      }
    }
  }
  
  par(mar = old.mar)
}
  
AlphaBFunc <- function(tau, b) {
  alpha <- 2 * (1 - b) / (tau * b - 2 * b + 2)
  return(alpha)
}


AlphaFunc <- function(tau, delta) {
  alpha <- 2 * delta / (1 - delta * tau + 2 * delta)
  alpha[alpha > 1] <- NA
  return(alpha)
}

PlotScnaCellFractions <- function(frac, seg.z, W, scna.ix, pal,
                                  bin.w, add=FALSE) {
  frac <- frac[scna.ix]
  seg.z <- seg.z[scna.ix]
  W <- W[scna.ix]
  
  colpal <- pal(1000)
  xlab <- "Fraction of cancer cells with alteration"
  x.max <- 1.0
  ix = ((!is.na(frac)) & (frac >= 0) & (frac < x.max))
   
  if (sum(ix) > 0) {
    PlotSeglenHist(frac[ix], W[ix], color.range=c(0,1), x.max=x.max,
                   bin.w=bin.w, order.by=seg.z[ix], use.pal=colpal,
                   data.whiskers=FALSE, xlab=xlab, x.axis=TRUE, add=add)
  } else {
    plot(0, type="n", bty="n", xlim=c(0,x.max), ylim=c(0,1), main="",
         xlab=xlab, ylab="", yaxt="n") 
  }
}
