
source('../../../scripts/funcs.R')
source('../../../scripts/libraries.R')
library(Seurat)

ov_sim_bulk = list()

##################### tcag #####################
OV_refPhi = readRDS('../../../../refPhi_cs/OV_refPhi.RDS')
ov_sim_bulk$TCGA = deconvBenchmarking::build_tcga_obj('OV',
                                                       to_use_genes = rownames(OV_refPhi@phi.cs))

#################### pseudobulk plus heter from Data_Izar2020_Ovarian ############
scExpr = Seurat::ReadMtx(mtx = '../../../../curated_dataset/OV/Data_Izar2020_Ovarian/10X/Exp_data_TPM_10X.mtx',
                         cells = '../../../../curated_dataset/OV/Data_Izar2020_Ovarian/10X/Cells_10X.csv',
                         features = '../../../../curated_dataset/OV/Data_Izar2020_Ovarian/10X/Genes_10X.txt',
                         feature.column = 1,skip.cell = 1,cell.sep = ',') # less memory usage compared with readMM()

scMeta = read.delim('../../../../curated_dataset/OV/Data_Izar2020_Ovarian/10X/Cells_10X.csv',sep = ',') %>% column_to_rownames('cell_name')
# this is a mixture of pre-treatment and post treatment samples

# remove cells with no assignment and Erythrocyte
to_keep = which(scMeta$cell_type!='' & scMeta$cell_type!='Erythrocyte')
scExpr = scExpr[,to_keep]
scMeta = scMeta[to_keep,]

ct_table = table(scMeta$cell_type)
simulated_frac = fracSimulator_Dirichlet(ct_table,n=30,dispersion_par = 0.0005,min.frac = 0.01)
scMeta$sample = paste0('sample_',scMeta$sample)

ov_sim_bulk$Izar2020 = create_pseudobulk_obj(scExpr,scMeta,unit = 'cpm',min.cells = 3,add_heter = T,simulated_frac = simulated_frac)

##################### pseudobulk plus heter from Data_Qian2020_Ovarian #############
scExpr = Seurat::ReadMtx(mtx = '../../../../curated_dataset/OV/Data_Qian2020_Ovarian/Exp_data_UMIcounts.mtx',
                         cells = '../../../../curated_dataset/OV/Data_Qian2020_Ovarian/Cells.csv',
                         features = '../../../../curated_dataset/OV/Data_Qian2020_Ovarian/Genes.txt',
                         feature.column = 1,skip.cell = 1,cell.sep = ',') # less memory usage compared with readMM()

scMeta = read.delim('../../../../curated_dataset/OV/Data_Qian2020_Ovarian/Cells.csv',sep = ',') %>% column_to_rownames('cell_name')
scMeta$sample = paste0('sample_',scMeta$sample)


to_keep = which(scMeta$cell_type!='')
scExpr = scExpr[,to_keep]
scMeta = scMeta[to_keep,]

ct_table = table(scMeta$cell_type)
simulated_frac = fracSimulator_Dirichlet(ct_table,n=20,dispersion_par = 0.0005,min.frac = 0.01)

ov_sim_bulk$Qian2020 = create_pseudobulk_obj(scExpr,scMeta,unit = 'UMI',min.cells = 3,
                                             add_heter = T,simulated_frac = simulated_frac)

saveRDS(ov_sim_bulk,file = 'sim_bulk.RDS')
