## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.

sample_trans_ccf_plot = function(alt_tab, mut_cols, grid, ccf_dens, sid, y_max) {
  x_max = length(grid)
  plot(0, type='n', bty='n', xaxt="n", main="", ylim=c(0, y_max), xlim=c(1, x_max), 
       xlab="Fraction of cancer cells with alteration", ylab="Density")
  axis(side=1, at=seq(0, x_max - 1, length=11), labels=seq(0, 1, by=0.1))
  title(main=sid, cex=par("cex"))
  
  trans = 1000 / nrow(alt_tab) 
  if (trans < 20) { 
    trans = 20 
  }
  if (trans > 255) { 
    trans = 255 
  }
  
  for(i in seq_len(nrow(ccf_dens))) {
    cr = col2rgb(mut_cols[i])
    use_col = rgb(cr[1, 1], cr[2, 1], cr[3, 1], trans, maxColorValue=255)
    
    xx = c(seq_len(ncol(ccf_dens)), rev(seq_len(ncol(ccf_dens))))
    yy = c(rep(0, ncol(ccf_dens)), rev(ccf_dens[i,]))
    polygon( xx, yy, border=FALSE, col=use_col, add=TRUE)
  }
}

GetMutBetaDensities <- function(mut.dat, n.grid=50) {
  mut.grid = matrix(0, nrow=nrow(mut.dat), ncol=n.grid)
  
  cov = mut.dat[, "alt"] + mut.dat[, "ref"]
  af = mut.dat[, "alt"] / cov
  
  grid.vals = seq_len(n.grid) / (n.grid + 1)
  
  for (i in seq_len(nrow(mut.grid))) {
    mut.grid[i, ] = dbeta(grid.vals, cov[i] * af[i] + 1, 
                          cov[i] * (1 - af[i]) + 1) * 1 / n.grid
  }
  
  return(mut.grid)
}

DrawMutBetaDensities <- function(beta.grid, pr.clonal, hz.del.flag, cols,
                                 draw.indv=TRUE, draw.total=TRUE) {
  n.grid <- ncol(beta.grid)
  grid.vals <- seq_len(n.grid) / (n.grid + 1)
  
  pal <- colorRampPalette(cols)
  colpal <- pal(1000)
  col.scale <- length(colpal)
  color.range <- c(0, 1)
  color.by <- pr.clonal
  pal.idx <- floor((color.by - color.range[1]) / (color.range[2] - color.range[1]) * 
    (col.scale - 1)) + 1
  mut.colors <- colpal[pal.idx]
  
  mut.colors[hz.del.flag] <- "navy"
  
  if (draw.indv) {
    for (i in seq_len(nrow(beta.grid))) {
      lines(grid.vals, beta.grid[i, ], col=mut.colors[i])
    }
  }
  
  ## Pr-weighted 
  clonal.grid <- matrix(0, nrow=nrow(beta.grid), ncol=ncol(beta.grid))
  sc.grid <- clonal.grid
  for (i in seq_len(nrow(beta.grid))) {
    clonal.grid[i, ] <- beta.grid[i, ] * pr.clonal[i]
    sc.grid[i, ] <- beta.grid[i, ] * (1 - pr.clonal[i])
  }
  
  if (draw.total) {
    lines(grid.vals, colSums(clonal.grid), col=cols[2], lty=4)
    lines(grid.vals, colSums(sc.grid), col=cols[1], lty=4)
  }
}

get_grid_combined_mut_densities = function(mut_pr, pr_clonal, grid, x_lim) {
  bin_w = x_lim / 100
  breaks=seq(0, x_lim, by=bin_w)
  mult_grid = breaks
  
  ## pr-weighted 
  clonal_grid = matrix(0, nrow=nrow(mut_pr), ncol=length(mult_grid))
  sc_grid = clonal_grid
  
  for (i in seq_len(nrow(mut_pr))) {
    x = grid[i, ] 
    y = mut_pr[i, ] * pr_clonal[i] 
    if (sum(!is.na(y)) > 2) {
      clonal_grid[i, ] = approx(x, y, xout=mult_grid)$y
    }
    y = mut_pr[i, ] * (1 - pr_clonal[i])
    
    if (sum(!is.na(y)) > 2) {
      sc_grid[i, ] = approx(x, y, xout=mult_grid)$y
    }
  }
  
  y_lim = max(c(colSums(clonal_grid, na.rm=TRUE), colSums(sc_grid, na.rm=TRUE)), na.rm=TRUE) 
  
  return(list(clonal_grid=clonal_grid, SC_grid=sc_grid, MULT_GRID=mult_grid, YLIM=y_lim))  
}


draw_grid_mut_densities = function(mut_pr, grid, pr_clonal, pr_cryptic_scna, x_lim, 
                                   xlab, cols, draw_indv=TRUE, draw_total=TRUE, 
                                   add=FALSE, y_lim=NA) {  
  pal = colorRampPalette(cols)
  colpal = pal(1000)
  
  col_scale = length(colpal)
  color_range = c(0, 1)
  color_by = pr_clonal
  pal_idx = floor((color_by - color_range[1]) / (color_range[2] - color_range[1]) * (col_scale-1)) + 1
  mut_colors = colpal[pal_idx]
  
  ix = pr_cryptic_scna > 0.5
  mut_colors[ix] = "mediumorchid2"
  
  # stack up variable-grid rescaled densities
  res = get_grid_combined_mut_densities(mut_pr, pr_clonal, grid, x_lim)
  clonal_grid = res$clonal_grid
  sc_grid = res$SC_grid 
  mult_grid= res$MULT_GRID
  
  if (!add) {
    if (is.na(y_lim)) { 
      y_lim = res$YLIM 
    }    
    plot( 0, type="n", bty="n", main="", xlab=xlab, ylab="Density", xlim=c(0,x_lim), 
          ylim=c(0, y_lim), las=1)
  }
  
  if (draw_indv) {
    for (i in seq_len(nrow(mut_pr))) {
      lines(grid[i, ], mut_pr[i, ], col=mut_colors[i])
    }
  }
  
  clonal_grid[is.na(clonal_grid)] = 0
  sc_grid[is.na(sc_grid)] = 0
  
  if (draw_total) {
    lines(mult_grid, colSums(clonal_grid), col=cols[2], lty=4)
    lines(mult_grid, colSums(sc_grid), col=cols[1], lty=4)
  }
}

PlotSomaticMutDensities <- function(mut.dat, seg.dat, mode.ix,
                                    mode.color, min.cov= 3, verbose=FALSE) {
  cov <- mut.dat[, "alt"] + mut.dat[, "ref"]
  ix <- cov > min.cov
  mut.dat <- mut.dat[ix, , drop=FALSE]
  cov <- cov[ix]
    
  if (nrow(mut.dat) == 0) {
    if (verbose) {
      print("No valid SNVs to plot")
    }
    for (i in 1:3) {
      frame()
    }
    return()
  }
  
  pr.clonal <- 1 - mut.dat[, "Pr_subclonal"]
  pr.cryptic.scna <- mut.dat[, "Pr_cryptic_SCNA"]

  draw.indv <- nrow(mut.dat) < 500
  n_grid = 300
  af_post_pr = GetMutBetaDensities(mut.dat, n_grid)
  grid = 1:n_grid / (n_grid + 1)
  grid_mat = matrix(grid, nrow=nrow(af_post_pr), ncol=ncol(af_post_pr), byrow=TRUE)
#  snv_cols = c("cyan4", "coral")
  snv_cols = c("maroon", "olivedrab4") 

  
  plot(0, type="n", bty="n", main="", xlab="SSNV allelic-fraction", ylab="Density", 
       xlim=c(0, 1), ylim=c(0, max(colSums(af_post_pr))), las=1 )
  hz.del.flag <- mut.dat[, "q_hat"] == 0
  DrawMutBetaDensities(af_post_pr, pr.clonal, hz.del.flag, snv_cols, draw.total=TRUE,
                       draw.indv=draw.indv)

  ## remove muts on HZdels
  mut.dat <- mut.dat[!hz.del.flag, , drop=FALSE]
  af_post_pr = af_post_pr[!hz.del.flag, , drop=FALSE]
  grid_mat = grid_mat[!hz.del.flag, , drop=FALSE]
  pr.clonal <- 1 - mut.dat[, "Pr_subclonal"]
  pr.cryptic.scna <- mut.dat[, "Pr_cryptic_SCNA"]
  
  if (nrow(mut.dat) == 0) {
    if (verbose) {
      print("No valid SNVs to plot")
    }
    for (i in 1:2) {
      frame()
    }
    return()
  }
  
  alpha <- mut.dat[1, "purity"]
  abline(v=alpha / 2, lty=2, col=mode.color)
  msg <- expression(hat(alpha) / 2)
  mtext(text=msg, side=3, at=alpha / 2, col=mode.color, line=0.2, cex=par("cex") * 
        par("cex.axis"))
      
  Q <- mut.dat[, "q_hat"]
  som.delta <- alpha / (2 * (1 - alpha) + alpha * Q)

  ## Plot rescaling to multiplicity
  mut.scale <- 1 / som.delta
  x = grid_mat * mut.scale
  y = af_post_pr * mut.scale^-1
  mult.xlim <- 2.5

  draw_grid_mut_densities(y, x, pr.clonal, pr.cryptic.scna, x_lim=mult.xlim,
                          cols=snv_cols, xlab="SSNV multiplicity", draw_total=TRUE,
                          draw_indv=draw.indv)
#  msg = "Rescaled to multiplicity"
#  mtext(text=msg, side=3, adj=0, cex=par("cex"))
  
  ## combined SCNAs and muts
  
  ## SNVs rescaling to cell-fraction
  ccf_grid = seq(0, 1, by=0.01)
  ccf_dens = matrix(NA, nrow=nrow(mut.dat), ncol=length(ccf_grid))
  for (i in seq_len(nrow(mut.dat))) {
    ccf_dens[i,] = calc_ccf_posterior_grid(mut.dat[i, "alt"], mut.dat[i, "ref"], 
                                           alpha, mut.dat[i, "q_hat"], ccf_grid)
  }
  
  ccf_grid_mat = matrix(ccf_grid, nrow=nrow(ccf_dens), ncol=ncol(ccf_dens), byrow=TRUE)
  
  ## Plot SCNA rescaling to cancer cell-fraction
  scna_cols = c("cyan4", "chocolate4") 
  sc_tab = seg.dat$mode.res$subclonal_SCNA_res$subclonal_SCNA_tab[mode.ix, , ]
  scna_ix = sc_tab[, "SCNA_ix"]
  scna_pr_subclonal = sc_tab[, "Pr_subclonal"]
  scna_pr_clonal = 1 - scna_pr_subclonal
  scna_pr_cryptic_scna = rep(0, length(pr.clonal))
  scna_ccf_hat = sc_tab[, "CCF_hat"]
  scna_ccf_dens = seg.dat$mode.res$subclonal_SCNA_res$CCF_dens[mode.ix, , ]
  scna_grid_mat = matrix(as.numeric(colnames(scna_ccf_dens)), nrow=nrow(scna_ccf_dens), 
                         ncol=ncol(scna_ccf_dens), byrow=TRUE)
  
  ix = !is.na(scna_ccf_hat) & scna_ix & sc_tab[, "subclonal_ix"]
  
  ## barplot
  if (sum(ix) > 0) {
    scna_res = get_grid_combined_mut_densities(scna_ccf_dens[ix, , drop=FALSE], 
                                               scna_pr_clonal[ix], 
                                               scna_grid_mat[ix, , drop=FALSE], x_lim=1)
    snv_res = get_grid_combined_mut_densities(ccf_dens, pr.clonal, ccf_grid_mat, x_lim=1)
    
    y_lim = max(scna_res$YLIM, snv_res$YLIM)

    draw_grid_mut_densities( scna_ccf_dens[ix, ,drop=FALSE], 
                             scna_grid_mat[ix, ,drop=FALSE], 
                             scna_pr_clonal[ix], 
                             scna_pr_cryptic_scna[ix], x_lim=1.0, cols=scna_cols, 
                             xlab="CCF", 
                             draw_total=TRUE, draw_indv=TRUE, add=FALSE, y_lim=y_lim )

    
#    alt_draw_grid_mut_densities(scna_ccf_dens[ix, ,drop=FALSE], 
#                                scna_grid_mat[ix, ,drop=FALSE], 
#                                scna_pr_clonal[ix], x_lim=1.0, 
#                                xlab="Fraction of cancer cells with alteration", 
#                                cols=c("black","black"), add=FALSE, y_lim=y_lim)
    
    add = TRUE
  } else{
    add = FALSE
  }

  draw_grid_mut_densities( ccf_dens, ccf_grid_mat, pr.clonal, scna_pr_cryptic_scna, x_lim=1.0, cols=snv_cols, xlab="CCF", draw_total=TRUE, draw_indv=draw.indv, add=add ) 

#  alt_draw_grid_mut_densities(ccf_dens, ccf_grid_mat, pr.clonal, x_lim=1.0, 
#                              cols=snv_cols, xlab="Fraction of cancer cells with alteration", 
#                              add=add)  

#  msg = "SSNVs and subclonal SCNAs"
#  mtext(text=msg, side=3, adj=0, cex=par("cex"))
   legend( x='topleft', legend=c("Clonal SSNVs", "Subclonal SSNVs", "Clonal SCNAs", "Subclonal SCNAs"), col=c( rev(snv_cols), rev(scna_cols) ), lty=1, lwd=1.5, bty="n" )
}


alt_draw_grid_mut_densities = function(mut_pr, grid, pr_clonal, x_lim, xlab, cols, 
                                       add=FALSE, y_lim=NA) { 
  ## add back clonal density
  m = ncol(mut_pr)
  
  mut_pr = mut_pr * (1 - pr_clonal)
  mut_pr[, m] = mut_pr[, m] + pr_clonal
  
  subclonal_ix = rep(NA, nrow(mut_pr))

  for (i in seq_len(nrow(mut_pr))) {
    sc_grid = grid[i, ] < 0.95
    subclonal_ix[i] = sum(mut_pr[i, sc_grid]) > 0.5
  }
  
  mut_colors = rep(NA, nrow(mut_pr))
  mut_colors[subclonal_ix] = cols[1]
  mut_colors[!subclonal_ix] = cols[2]
  
  if (!add) {
    if (is.na(y_lim)) { 
      y_lim = max(mut_pr) 
    }
    
    plot(0, type="n", bty="n", main="", xlab=xlab, ylab="Density", xlim=c(0, x_lim), 
         ylim=c(0, y_lim), las=1)
  }
  
  for (i in seq_len(nrow(mut_pr))) {
    lines( grid[i,], mut_pr[i,], col=mut_colors[i] )
  }
}




