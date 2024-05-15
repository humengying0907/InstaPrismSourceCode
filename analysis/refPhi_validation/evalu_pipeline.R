#!/usr/bin/env Rscript
library(argparse)
parser <- ArgumentParser(description = 'InstaPrism reference evaluation pipeline',add_help = T)

parser$add_argument("--tumorType",'-t',type="character",metavar="<character>",help = "Folder name of the bulk dataset for evaluation. The 'sim_bulk.RDS' file located within '../bulk_simulation/<tumorType>' will be utilized for deconvolution." )
parser$add_argument("--refName",'-n',type="character", nargs='+',metavar="<character>",help = "Name of the reference to test. The reference stored in '../../refPhi/<refName>_refPhi.RDS' will be loaded as the reference. To test multiple references, separate each name by a blank space.")
parser$add_argument("--output",'-o',type="character",required = F,metavar="<character>",help = "Output folder name. A 'performance/<output>/' directory will be created to store the deconvolution results and plots. If not specified,  a 'performance/<tumorType>/' directory will be created to store the results" )
parser$add_argument("--niter", type="integer",required = F,default = 400,metavar="<integer>",help = 'Number of iterations for InstaPrism() function. [default: 400]')
parser$add_argument("--ncore", type="integer",required = F,default = 16, metavar="<integer>",help = 'Number of threads. [default: 16]')
parser$add_argument('--updateReference',type="logical",required = F,default = F, metavar="<logical>",help = 'A logical variable to determine whether to include updated reference in the evaulation. [default: FALSE]')
parser$add_argument("--key", type="character",required = F,metavar="<character>",help = 'Name of the malignant cell type in the reference, required only when updateReference = TRUE. When evaluating multiple references, these references need to have the same key.') 
parser$add_argument('--saveDeconvRes',type="logical",required = F,default = F,metavar="<logical>", help = 'A logical variable to determine whether to save the deconvolution results. [default: FALSE]')

args <- parser$parse_args()

tumorType = args$tumorType
key = args$key
niter = args$niter
ncore = args$ncore
updateReference = args$updateReference
saveDeconvRes = args$saveDeconvRes

if(is.null(args$output)){
  output = tumorType
}else{
  output = args$output
}

if(updateReference){
  if(is.null(key)){
    stop('need to provide a key when updateReference = TRUE')
  }
}

source('../../scripts/funcs.R')
source('../../scripts/libraries.R')

# load sim_bulk
sim_bulk = readRDS(paste0('../bulk_simulation/',tumorType,'/sim_bulk.RDS'))


refPhi_list = list()
for(refName in args$refName){
  refPhi_list[[refName]] = readRDS(paste0('../../refPhi/',refName,'_refPhi.RDS'))
}

# InstaPrism deconvolution
deconv_res = multi_deconv(sim_bulk,refPhi_list,key,updateReference,niter,ncore)

# create a folder to store performance data associated with this reference
if(!dir.exists(paste0('performance/',output))){
  dir.create(paste0('performance/',output))
}

if(saveDeconvRes){
  saveRDS(deconv_res,file = paste0('performance/',output,'/deconv_res.RDS'))
}

# theta performance
theta_res = extract_theta(deconv_res)
saveRDS(theta_res,file = paste0('performance/',output,'/theta_res.RDS'))

theta_performance = list()
for(bulk in names(theta_res)){
  theta_performance[[bulk]] = theta_evalu(theta_res[[bulk]],sim_bulk[[bulk]]$simulated_frac)
}
saveRDS(theta_performance,file = paste0('performance/',output,'/theta_performance.RDS'))


# visualize theta performance
for(bulk in names(sim_bulk)){
  print(bulk)
  n_ct = nrow(theta_performance[[bulk]][[1]][["summ"]])
  n_ref = length(theta_performance[[bulk]])
  
  theta_plot_export(theta_performance[[bulk]],fig_width = n_ct * 2.5, fig_height = n_ref * 2.5,
                    plot_file_name = paste0('performance/',output,'/',bulk,'_theta_performance.pdf'))
}


