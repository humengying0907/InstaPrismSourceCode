
source('../../../scripts/funcs.R')
source('../../../scripts//libraries.R')
library(Seurat)

brca_sim_bulk = list()

##################### tcag #####################
BRCA_refPhi = readRDS('../../../refPhi/BRCA_refPhi.RDS')
brca_sim_bulk$TCGA = deconvBenchmarking::build_tcga_obj('BRCA',
                                                         to_use_genes = rownames(BRCA_refPhi@phi.cs))


############ pseudobulk plus heter simulated bulk samples from Qian2020 ##########
scExpr = Seurat::ReadMtx(mtx = '../../../curated_dataset/BRCA/Data_Qian2020_Breast/Exp_data_UMIcounts.mtx',
                         cells = '../../../curated_dataset/BRCA/Data_Qian2020_Breast/Cells.csv',
                         features = '../../../curated_dataset/BRCA/Data_Qian2020_Breast/Genes.txt',
                         feature.column = 1,skip.cell = 1,cell.sep = ',') # less memory usage compared with readMM()
scMeta = read.delim('../../../curated_dataset/BRCA/Data_Qian2020_Breast/Cells.csv',sep =',')
scMeta = scMeta %>% column_to_rownames('cell_name')

simulated_frac = fracSimulator_Dirichlet(table(scMeta$cell_type),n=30,dispersion_par = 0.0005,min.frac = 0.01)
brca_sim_bulk$Qian2020 = create_pseudobulk_obj(scExpr,scMeta,unit = 'UMI',add_heter = T,simulated_frac = simulated_frac)


saveRDS(brca_sim_bulk,file = 'sim_bulk.RDS')



