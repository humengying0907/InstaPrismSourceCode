
source('../../../scripts/funcs.R')
source('../../../scripts/libraries.R')
library(Seurat)

gbm_sim_bulk = list()

##################### tcag #####################
GBM_refPhi = readRDS('../../../refPhi/GBM_refPhi.RDS')
gbm_sim_bulk$TCGA = deconvBenchmarking::build_tcga_obj('GBM',
                                                       to_use_genes = rownames(GBM_refPhi@phi.cs))
####### pseudobulk from Data_Neftel2019_Brain ##########
scExpr = Seurat::ReadMtx(mtx = '../../../../curated_dataset/GBM/Data_Neftel2019_Brain/SmartSeq2/exp_data_TPM.mtx',
                         cells = '../../../../curated_dataset/GBM/Data_Neftel2019_Brain/SmartSeq2/Cells.csv',
                         features = '../../../../curated_dataset/GBM/Data_Neftel2019_Brain/SmartSeq2/genes.txt',
                         feature.column = 1,skip.cell = 1,cell.sep = ',') # less memory usage compared with readMM()
scMeta = read.delim('../../../../curated_dataset/GBM/Data_Neftel2019_Brain/SmartSeq2/Cells.csv',sep = ',') %>% column_to_rownames('cell_name')
gbm_sim_bulk$Neftel2019 = create_pseudobulk_obj(scExpr,scMeta,unit = 'cpm')

#################### pseudobulk from Yuan2018 ##################
scExpr = Seurat::ReadMtx(mtx = '../../../../curated_dataset/GBM/Data_Yuan2018_Brain/Exp_data_UMIcounts.mtx',
                         cells = '../../../../curated_dataset/GBM/Data_Yuan2018_Brain/Cells.csv',
                         features = '../../../../curated_dataset/GBM/Data_Yuan2018_Brain/Genes.txt',
                         feature.column = 1,skip.cell = 1,cell.sep = ',') # less memory usage compared with readMM()
scMeta = read.delim('../../../../curated_dataset/GBM/Data_Yuan2018_Brain/Cells.csv',sep = ',') %>% column_to_rownames('cell_name')

ct_table = table(scMeta$cell_type)

simulated_frac = fracSimulator_Dirichlet(ct_table,n=50,dispersion_par = 0.0003,min.frac = 0.01)

gbm_sim_bulk$Yuan2018 = create_pseudobulk_obj(scExpr,scMeta,unit = 'UMI',min.cells = 3,add_heter = T,simulated_frac = simulated_frac)
saveRDS(gbm_sim_bulk,file = 'sim_bulk.RDS')
