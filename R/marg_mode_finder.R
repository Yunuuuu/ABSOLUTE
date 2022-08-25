## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.


GetAlphaAndTau = function(b, delta) {
   alpha = (2 * delta) / (2 * delta + b)
   tau = -(b - 1) / delta

   return(list(alpha=alpha, tau=tau)) 
}

FitSample = function(seg.obj, mut, Q, pi_theta_qz, sigma.h, tau.dom,
                      sigma.h.dom, chr.arms.dat, verbose=FALSE) {
  kThetaQz = rep(1 / (Q + 1), Q + 1)

  obs = seg.obj[["obs.scna"]]
  
  res = FindLocationModes(obs, mut, Q, kThetaQz, sigma.h,
                           tau.dom, verbose=verbose)
  
  if (!is.null(res[["mode.flag"]])) { 
      return(list("mode.flag"=res[["mode.flag"]]))
   } else {
     mode.tab = res[["mode.tab"]]
   }

   if (is.na(mode.tab))  {
      return(list(mode.flag="ERROR"))
   }

   n.modes = nrow(mode.tab)
   theta.qz.hat = matrix(NA, nrow=n.modes, ncol=Q+1)
   theta.q.tab = array(NA, dim=c(n.modes, Q))

  if (verbose) {
    print(paste("Optimizing LL(data, theta.qz, sigma.h | comb) for ",
                n.modes, " modes: ", sep=""))
  }
  
  seg.z.tab = array(NA, dim=c(n.modes, length(obs[["W"]])))
  seg.qz.tab = array(NA, dim=c(n.modes, length(obs[["W"]]), Q+1))
  seg.q.tab = array(NA, dim=c(n.modes, length(obs[["W"]]), Q))
  
  if (obs[["data.type"]] == "ALLELIC") {
    ab.tab = array(NA, dim=c(n.modes, Q+1))
    chr.arm.tab = array(NA, dim=c(n.modes, 2, nrow(chr.arms.dat), Q))
  }
  if (obs[["data.type"]] == "TOTAL") {
    ab.tab = NULL
    chr.arm.tab = array(NA, dim=c(n.modes, 1, nrow(chr.arms.dat), Q))
  }
  
  dimnames(chr.arm.tab)[[3]] = rownames(chr.arms.dat)
  mode.tab = cbind(mode.tab, NA, NA, NA, NA, NA, NA, NA, NA)
  colnames(mode.tab)[c((ncol(mode.tab)-7) : ncol(mode.tab))] = c("genome mass", "sigma.h.hat", "theta.z.hat", "frac.het", "log.partition", "entropy", "clust_LL", "post_LL")
  
  mode.tab = cbind(mode.tab, "somatic.mut.ll"=NA) 

  for (i in 1:n.modes) {
    delta = mode.tab[i, "delta"]
    b = mode.tab[i, "b"]
    ## can be NA for non-array platforms
    obs[["error.model"]][["fit.at"]] = mode.tab[i, "AT"]  
    comb =  GetCopyRatioComb(Q, delta, b, obs[["error.model"]])

    ## optimize
    res = OptThetaQzSigmaH(obs, comb, sigma.h.dom, pi_theta_qz, verbose=verbose)

    mode.tab[i, "sigma.h.hat"] = res[["sigma.h.hat"]]
    mode.tab[i, "log.partition"] = res[["LL"]]
    ## initial version with data-fit term only
    mode.tab[i, "post_LL"] = res[["LL"]]  
    lambda.qz.res = res[["lambda.qz.res"]]
    mode.tab[i, "theta.z.hat"] = res$theta_qz_hat[(Q + 1)]
    
    ## posterior distribution over segment CNs
    res = GetSegQzPost(obs[["d.tx"]], obs[["d.stderr"]],
                        obs[["W"]], comb, NA, mode.tab[i,"sigma.h.hat"],
                        theta.qz=res[["theta_qz_hat"]])
    seg.z.tab[i, ] = res[["QZ"]][, (Q+1)]
    seg.qz.tab[i, , ] = res[["QZ"]]
    seg.q.tab[i, , ] = res[["Q"]]
    theta.q.tab[i, ] = colSums(res[["Q"]] * obs[["W"]])

    mode.tab[i, "theta.z.hat"] = sum(seg.qz.tab[i, , (Q + 1)] * obs$W)
    
    ## compute % non-clonal genome
    mode.tab[i, "frac.het"] = sum(obs[["W"]] * seg.z.tab[i,])  

    ## calculate allelic-balance
    if (obs[["data.type"]] == "ALLELIC") {
      ab.tab[i, ] = CalcAbDistr(obs, res[["QZ"]])
      mode.tab[i,"genome mass"] = 2 * sum(c((1:Q)-1) * colSums(res[["Q"]] * obs[["W"]]))
    }
    if (obs[["data.type"]] == "TOTAL") {
      mode.tab[i,"genome mass"] = 1 * sum(c((1:Q)-1) * colSums(res[["Q"]] * obs[["W"]]))
    }

    ## weighted entropy average over segs
    mode.tab[i, "entropy"] = CalcFitEntropy(obs, res[["QZ"]])
    
    chr.arm.tab[i, , , ] = CalcChrArmDistr(seg.obj, res[["Q"]], chr.arms.dat)
  }

  return(list(mode.tab=mode.tab, mode.posts=NA,
              theta.q.tab=theta.q.tab, b=theta.qz.hat,
              seg.z.tab=seg.z.tab, seg.qz.tab=seg.qz.tab,
              seg.q.tab=seg.q.tab, ab.tab=ab.tab, chr.arm.tab=chr.arm.tab,
              mode.flag=NA))
}

WeighSampleModes = function(mode.res) {
  mode.ll = mode.res[["mode.tab"]][, "post_LL"]
  mx = max(mode.ll)
  nf = sum(exp(mode.ll - mx))
  dens = exp(mode.ll - mx) / nf
  mode.res[["mode.tab"]] = cbind(mode.res[["mode.tab"]], dens)
  
  if (all(is.finite(dens))) {
    ix = sort(mode.res[["mode.tab"]][,"dens"], decreasing=TRUE,
               index.return=TRUE)[["ix"]]
  } else {
    ix = c(1:nrow(mode.res[["mode.tab"]]))
  }
  
  mode.res = ReorderModeRes(mode.res, ix)
  
  return(mode.res)
}

ReorderModeRes = function(mode.res, ix, DROP=FALSE) {
   mode.res[["mode.tab"]] = mode.res[["mode.tab"]][ix,, drop=DROP]
   mode.res[["seg.z.tab"]] = mode.res[["seg.z.tab"]][ix,, drop=DROP]
   mode.res[["seg.qz.tab"]] = mode.res[["seg.qz.tab"]][ix,,, drop=DROP]
   mode.res[["seg.q.tab"]] = mode.res[["seg.q.tab"]][ix,,, drop=DROP]

   ## only exists for allelic data
   if (!is.null(mode.res[["ab.tab"]]))  {
      mode.res[["ab.tab"]] = mode.res[["ab.tab"]][ix, , drop=DROP]
   }

   mode.res[["theta.q.tab"]] = mode.res[["theta.q.tab"]][ix, ,drop=DROP]
   mode.res[["theta.qz.hat"]] = mode.res[["theta.qz.hat"]][ix, ,drop=DROP]
   mode.res[["chr.arm.tab"]] = mode.res[["chr.arm.tab"]][ix ,,, , drop=DROP]
   mode.res[["mode.clust.p"]] = mode.res[["mode.clust.p"]][ix , , drop=DROP]

   mode.res$subclonal_SCNA_res$subclonal_SCNA_tab = mode.res$subclonal_SCNA_res$subclonal_SCNA_tab[ix, , , drop=DROP]
   mode.res$subclonal_SCNA_res$CCF_dens = mode.res$subclonal_SCNA_res$CCF_dens[ix, , , drop=DROP]
   
   ## only exists if MAF supplied.
   if (!is.null(mode.res[["muts.post.prs"]])) {
      mode.res[["muts.post.prs"]] = mode.res[["muts.post.prs"]][,,ix, drop=DROP]
   }
   mode.res[["mode.posts"]] = mode.res[["mode.posts"]][ix]

   return(mode.res)
}

## dont call non-aneuploid if mut_dat is present
GetCallStatus = function(mode.res, seg.w) {
  status = "called"
  q.tab = mode.res[["seg.qz.tab"]][1, , ]
  q.tab = q.tab[, c(1:  (ncol(q.tab)) )]
  
  Q = ncol(q.tab) - 1
  max.q = (apply(q.tab, 1, which.max))

  peak_masses = rep(0, Q)
  for (i in 1:Q) {
    ix = which(max.q==i) 
    peak_masses[i] = sum(seg.w[ix])
  }
  
  b = mode.res[["mode.tab"]][1,"b"]
  
  ## don't count 0-CN state if b is too small - could be due to IBD, not LOH
  if (peak_masses[1] < 0.01 & b < 0.15) {
    peak_masses[1] = 0
  }

  six = order(peak_masses, decreasing=TRUE)

  if (peak_masses[six[3]] < 0.0001) {
    if (peak_masses[six[2]] < 0.0001) { 
      if (mode.res[["mode.tab"]][1,"sigma.h.hat"] < 0.02 &
          !is.null(mode.res[["muts.post.prs"]])) { 
        status = "non-aneuploid" 
      } else {
        status = "low purity"
      }
    }
  }

   if (mode.res[["mode.tab"]][1, "entropy"] > 0.2) {
     status = "high entropy"
   }

   if (mode.res[["mode.tab"]][1, "frac.het"] > 0.2) {
     status = "high non-clonal"
   }

   return(status)
}
