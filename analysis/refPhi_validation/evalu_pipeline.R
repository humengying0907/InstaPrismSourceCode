#!/usr/bin/env Rscript
library(argparse)
parser <- ArgumentParser(add_help = F)

parser$add_argument("--tumorType",'-t',type="character")
parser$add_argument("--key", type="character") # set to None if no key applied
parser$add_argument("--niter", type="integer",required = F,default = 300)
parser$add_argument("--deconv_TCGA", type="logical",required = F,default = T)

args <- parser$parse_args()

tumorType = args$tumorType
key = args$key
niter = args$niter
deconv_TCGA = args$deconv_TCGA

source('../funcs.R')
source('../libraries.R')
source('../libraries_plot.R')

# load refPhi_cs object
refPhi_cs = readRDS(paste0('../../refPhi_cs/',tumorType,'_refPhi.RDS'))

# load sim_bulk
sim_bulk = readRDS(paste0('bulk/',tumorType,'/sim_bulk.RDS'))

if(!deconv_TCGA){
  sim_bulk = sim_bulk[names(sim_bulk)!= 'TCGA']
}


if(key == 'None'){
  ct.update.list = list()
  ct.update.list[[1]] =  'all'
  key = NA
}else{
  ct.update.list = list() # specify a list of cell.types.to.update parameter used for InstaPrism_update() function
  ct.update.list[[1]] =  NULL # update mal only
  ct.update.list[[2]] =  'all'
}


# InstaPrism deconvolution
deconv_res = multi_deconv(sim_bulk,refPhi_cs,key,ct.update.list,niter)
if(!dir.exists(paste0('performance/',tumorType))){
  dir.create(paste0('performance/',tumorType))
}

saveRDS(deconv_res,file = paste0('performance/',tumorType,'/deconv_res.RDS'))

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

# Rscript evalu_pipeline.R -t LUAD --key Malignant --niter 400
# Rscript evalu_pipeline.R -t GBM --key Malignant --niter 400 
# Rscript evalu_pipeline.R -t OV --key malignant --niter 400 --deconv_TCGA FALSE 
# Rscript evalu_pipeline.R -t RCC --key Malignant --niter 400 --deconv_TCGA FALSE 

