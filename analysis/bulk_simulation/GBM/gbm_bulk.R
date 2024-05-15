
source('../../../scripts/funcs.R')
source('../../../scripts/libraries.R')
library(Seurat)

gbm_sim_bulk = list()

##################### tcag #####################
GBM_refPhi = readRDS('../../../refPhi/GBM_refPhi.RDS')
gbm_sim_bulk$TCGA = deconvBenchmarking::build_tcga_obj('GBM',
                                                       to_use_genes = rownames(GBM_refPhi@phi.cs))
####### pseudobulk from Data_Neftel2019_Brain ##########
dataset_path = '../../../../../curated_dataset/GBM/Data_Neftel2019_Brain/SmartSeq2/'


scExpr = Seurat::ReadMtx(mtx = paste0(dataset_path,'exp_data_TPM.mtx'),
                         cells = paste0(dataset_path,'Cells.csv'),
                         features = paste0(dataset_path,'genes.txt'),
                         feature.column = 1,skip.cell = 1,cell.sep = ',') 

scMeta = read.delim(paste0(dataset_path,'Cells.csv'),sep = ',') %>% column_to_rownames('cell_name')



n_pseudo = length(unique(scMeta$sample))
n_simu = 50-n_pseudo

simulated_frac = fracSimulator_Dirichlet(table(scMeta$cell_type),n=n_simu,dispersion_par = 0.0005,min.frac = 0.01)
gbm_sim_bulk$Neftel2019 = create_pseudobulk_obj(scExpr,scMeta,unit = 'cpm',add_heter = T,simulated_frac = simulated_frac)


#################### pseudobulk from Yuan2018 ##################
dataset_path = '../../../../../curated_dataset/GBM/Data_Yuan2018_Brain/'

scExpr = Seurat::ReadMtx(mtx = paste0(dataset_path,'Exp_data_UMIcounts.mtx'),
                         cells = paste0(dataset_path,'Cells.csv'),
                         features = paste0(dataset_path,'Genes.txt'),
                         feature.column = 1,skip.cell = 1,cell.sep = ',') 
scMeta = read.delim(paste0(dataset_path,'Cells.csv'),sep = ',') %>% column_to_rownames('cell_name')

ct_table = table(scMeta$cell_type)
n_pseudo = length(unique(scMeta$sample))
n_simu = 50-n_pseudo

# increase prior for 'Oligodendrocyte','T_cell' to avoid too small fraction in the simulation
ct_table[c('Oligodendrocyte','T_cell','Pericyte','Endothelial')] = ct_table[c('Oligodendrocyte','T_cell','Pericyte','Endothelial')]+500

simulated_frac = fracSimulator_Dirichlet(ct_table,n=n_simu,dispersion_par = 0.0003,min.frac = 0.01)

gbm_sim_bulk$Yuan2018 = create_pseudobulk_obj(scExpr,scMeta,unit = 'UMI',min.cells = 3,
                                              add_heter = T,simulated_frac = simulated_frac,
                                              min_chunkSize = 50)
saveRDS(gbm_sim_bulk,file = 'sim_bulk.RDS')
