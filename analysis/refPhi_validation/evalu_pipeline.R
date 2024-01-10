#!/usr/bin/env Rscript
library(argparse)
parser <- ArgumentParser(add_help = F)

parser$add_argument("--tumorType",'-t',type="character")
parser$add_argument("--key", type="character") 
parser$add_argument("--niter", type="integer",required = F,default = 400)
parser$add_argument("--ncore", type="integer",required = F,default = 1)
parser$add_argument('--updateReference',type="logical",required = F,default = T)
parser$add_argument('--saveDeconvRes',type="logical",required = F,default = F)

args <- parser$parse_args()

tumorType = args$tumorType
key = args$key
niter = args$niter
ncore = args$ncore
updateReference = args$updateReference
saveDeconvRes = args$saveDeconvRes

source('../../scripts/funcs.R')
source('../../scripts/libraries.R')


# load refPhi_cs object
refPhi_cs = readRDS(paste0('../../refPhi/',tumorType,'_refPhi.RDS'))

# load sim_bulk
sim_bulk = readRDS(paste0('../bulk_simulation/',tumorType,'/sim_bulk.RDS'))

# InstaPrism deconvolution
deconv_res = multi_deconv(sim_bulk,refPhi_cs,key,updateReference,niter,ncore)
if(!dir.exists(paste0('performance/',tumorType))){
  dir.create(paste0('performance/',tumorType))
}

if(saveDeconvRes){
  saveRDS(deconv_res,file = paste0('performance/',tumorType,'/deconv_res.RDS'))
}

# theta performance
theta_res = extract_theta(deconv_res)
theta_performance = theta_evalu(theta_res,sim_bulk)
saveRDS(theta_performance,file = paste0('performance/',tumorType,'/theta_performance.RDS'))


# visualize theta performance
for(name in names(sim_bulk)){
  print(name)
  n_ct = nrow(theta_performance[[name]][["deconvEvalu"]][["initial"]][["summ"]])
  theta_plot_export(theta_performance[[name]],fig_width = n_ct * 2.5, fig_height = 8,plot_file_name = paste0('performance/',tumorType,'/',name,'_theta_performance.pdf'))
}


