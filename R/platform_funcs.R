## The Broad Institute
## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2012) by the
## Broad Institute/Massachusetts Institute of Technology. All rights are
## reserved.
##
## This software is supplied without any warranty or guaranteed support
## whatsoever. Neither the Broad Institute nor MIT can be responsible for its
## use, misuse, or functionality.


## FIXME: Some funcs don't seem to be used at all:
##   - GetMutsScna
##   - CalcSampleMutsPostPr
##   - get_cr_grid_from_ccf_grid
set_allelic_funcs = function() {
  MakeSegObj <<- AllelicMakeSegObj
  ExtractSampleObs <<- AllelicExtractSampleObs
  GetMutSegIx <<- AllelicGetMutSegIx
  get_muts_nearest_clonal_scna <<- allelic_get_muts_nearest_clonal_scna
  CalcSampleMutsPostPr <<- AllelicCalcSampleMutsPostPr
  Atten <<- AllelicAtten
  TxData <<- AllelicTxData
  InvTxData <<- AllelicInvTxData
  InvAtten <<- AllelicInvAtten
  GetScnaStderrGridDensity <<- AllelicGetScnaStderrGridDensity
  CalcChrArmDistr <<- AllelicCalcChrArmDistr
  GetCopyRatioComb <<- AllelicGetCopyRatioComb
  get_cr_grid_from_ccf_grid <<- allelic_get_cr_grid_from_ccf_grid
  GetAbsSegDat <<- AllelicGetAbsSegDat
  calc_sample_muts_on_subclonal_scna <<- allelic_calc_sample_muts_on_subclonal_scna
  get_subclonal_scna_mut_ix <<- allelic_get_subclonal_scna_mut_ix  
} 

set_total_funcs = function() {
  MakeSegObj <<- total_make_seg_obj
  ExtractSampleObs <<- total_extract_sample_obs
  GetMutSegIx <<- total_get_mut_seg_ix    
  get_muts_nearest_clonal_scna<<- total_get_muts_nearest_clonal_scna
  CalcSampleMutsPostPr <<- total_calc_sample_muts_post_pr
  Atten <<- total_atten
  TxData <<- total_tx_data
  InvTxData <<- total_inv_tx_data
  InvAtten <<- total_inv_atten
  GetScnaStderrGridDensity <<- total_get_scna_stderr_grid_density
  CalcChrArmDistr <<- total_calc_chr_arm_distr
  GetCopyRatioComb <<- total_get_copy_ratio_comb
  get_cr_grid_from_ccf_grid <<- total_get_cr_grid_from_ccf_grid
  GetAbsSegDat <<- total_get_abs_seg_dat
  calc_sample_muts_on_subclonal_scna <<- total_calc_sample_muts_on_subclonal_scna  
  get_subclonal_scna_mut_ix <<- total_get_subclonal_scna_mut_ix
}