
source('../../../scripts/funcs.R')
source('../../../scripts/libraries.R')
library(Seurat)

ov_sim_bulk = list()

##################### tcag #####################
OV_refPhi = readRDS('../../../../refPhi_cs/OV_refPhi.RDS')
ov_sim_bulk$TCGA = deconvBenchmarking::build_tcga_obj('OV',
                                                       to_use_genes = rownames(OV_refPhi@phi.cs))

#################### pseudobulk plus heter from Data_Izar2020_Ovarian ############
dataset_path = '../../../../../curated_dataset/OV/Data_Izar2020_Ovarian/10X/'

scExpr = Seurat::ReadMtx(mtx = paste0(dataset_path,'Exp_data_TPM_10X.mtx'),
                         cells = paste0(dataset_path,'Cells_10X.csv'),
                         features = paste0(dataset_path,'Genes_10X.txt'),
                         feature.column = 1,skip.cell = 1,cell.sep = ',') 

scMeta = read.delim(paste0(dataset_path,'Cells_10X.csv'),sep = ',') %>% column_to_rownames('cell_name')
# this is a mixture of pre-treatment and post treatment samples

# remove cells with no assignment and Erythrocyte since it only contains few cells
to_keep = which(scMeta$cell_type!='' & scMeta$cell_type!='Erythrocyte')
scExpr = scExpr[,to_keep]
scMeta = scMeta[to_keep,]
scMeta$sample = paste0('sample_',scMeta$sample)

ct_table = table(scMeta$cell_type)
n_pseudo = length(unique(scMeta$sample))
n_simu = 50-n_pseudo

ct_table[c(1,2,5,6)] = ct_table[c(1,2,5,6)]*3

# increase prior to avoid too small fraction in the simulation
simulated_frac = fracSimulator_Dirichlet(ct_table,n=n_simu,dispersion_par = 0.0005,min.frac = 0.01)
ov_sim_bulk$Izar2020 = create_pseudobulk_obj(scExpr,scMeta,unit = 'UMI',min.cells = 3,
                                             add_heter = T,simulated_frac = simulated_frac,
                                             min_chunkSize = 50)

##################### pseudobulk plus heter from Data_Qian2020_Ovarian #############
dataset_path = '../../../../../curated_dataset/OV/Data_Qian2020_Ovarian/'

scExpr = Seurat::ReadMtx(mtx = paste0(dataset_path,'Exp_data_UMIcounts.mtx'),
                         cells = paste0(dataset_path,'Cells.csv'),
                         features = paste0(dataset_path,'Genes.txt'),
                         feature.column = 1,skip.cell = 1,cell.sep = ',') 

scMeta = read.delim(paste0(dataset_path,'Cells.csv'),sep = ',') %>% column_to_rownames('cell_name')

scMeta$sample = paste0('sample_',scMeta$sample)


to_keep = which(scMeta$cell_type!='')
scExpr = scExpr[,to_keep]
scMeta = scMeta[to_keep,]

ct_table = table(scMeta$cell_type)

n_pseudo = length(unique(scMeta$sample))
n_simu = 50-n_pseudo

simulated_frac = fracSimulator_Dirichlet(ct_table,n=n_simu,dispersion_par = 0.0005,min.frac = 0.01)

ov_sim_bulk$Qian2020 = create_pseudobulk_obj(scExpr,scMeta,unit = 'UMI',min.cells = 3,
                                             add_heter = T,simulated_frac = simulated_frac,min_chunkSize = 50)

saveRDS(ov_sim_bulk,file = 'sim_bulk.RDS')
