
source('../../../scripts/funcs.R')
source('../../../scripts/libraries.R')
library(Seurat)

RCC_refPhi = readRDS('../../../refPhi/RCC_refPhi.RDS')
rcc_sim_bulk = list()

########################### TCGA ############################
rcc_sim_bulk$TCGA = deconvBenchmarking:::build_tcga_obj('KIRC',
                                                        to_use_genes = rownames(RCC_refPhi@phi.cs),
                                                        purity_methods =  c('ESTIMATE', 'ABSOLUTE', 'LUMP', 'IHC', 'CPE','ABSOLUTE_GDC'))

################## pseudobulk from Zhang2021_Kidney ###################
scExpr = Seurat::ReadMtx(mtx = '../../../../curated_dataset/RCC/Data_Zhang2021_Kidney/Exp_data_UMIcounts.mtx',
                         cells = '../../../../curated_dataset/RCC/Data_Zhang2021_Kidney/Cells.csv',
                         features = '../../../../curated_dataset/RCC/Data_Zhang2021_Kidney/Genes.txt',
                         feature.column = 1,skip.cell = 1,cell.sep = ',') # less memory usage compared with readMM()

scMeta = read.delim('../../../../curated_dataset/RCC/Data_Zhang2021_Kidney/Cells.csv',sep = ',') %>% column_to_rownames('cell_name')

# use ccRCC only
keep_id = which(scMeta$cell_type!='' & scMeta$disease == 'Clear_Cell_RCC')
scExpr = scExpr[,keep_id]
scMeta = scMeta[keep_id,]


ct_table = table(scMeta$cell_type)
simulated_frac = fracSimulator_Dirichlet(ct_table,n=50,dispersion_par = 0.0003,min.frac = 0.01)
rcc_sim_bulk$Zhang2021 = create_pseudobulk_obj(scExpr,scMeta,
                                               unit = 'UMI',
                                               min.cells = 3,
                                               add_heter = T,
                                               simulated_frac = simulated_frac,
                                               min_chunkSize = 10)

#### pseuodbulk from Young2018_Kidney #######
scExpr = Seurat::ReadMtx(mtx = '../../../../curated_dataset/RCC/Data_Young2018_Kidney/Exp_data_UMIcounts.mtx',
                         cells = '../../../../curated_dataset/RCC/Data_Young2018_Kidney/Cells.csv',
                         features = '../../../../curated_dataset/RCC/Data_Young2018_Kidney/Genes.txt',
                         feature.column = 1,skip.cell = 1,cell.sep = ',') # less memory usage compared with readMM()

scMeta = read.delim('../../../../curated_dataset/RCC/Data_Young2018_Kidney/Cells.csv',sep = ',') %>% column_to_rownames('cell_name')
scMeta$disease = gsub("^(.*?)_.*$", "\\1", scMeta$sample)


to_keep = which(scMeta$cell_QCpass == T & scMeta$cell_type!='' 
                & scMeta$cell_type!= 'Erythroblast' 
                & scMeta$cell_type!= 'Neutrophil'
                & scMeta$disease %in% c('RCC1','RCC2','RCC3'))

scExpr = scExpr[,to_keep]
scMeta = scMeta[to_keep,]

# before generating pseudobulk samples, further remove samples with few cells
to_use = which(!scMeta$sample %in% c('RCC3_Kid_T_ldc_1_3','RCC3_Kid_T_ldc_1_4'))
scExpr = scExpr[,to_use]
scMeta = scMeta[to_use,]

ct_table = table(scMeta$cell_type)
simulated_frac = fracSimulator_Dirichlet(ct_table,n=50,dispersion_par = 0.0003,min.frac = 0.001)

rcc_sim_bulk$Young2018 = create_pseudobulk_obj(scExpr,scMeta,
                                               unit = 'UMI',
                                               min.cells = 10,
                                               add_heter = T,
                                               simulated_frac = simulated_frac,
                                               # increase min_chunksize to avoid sparsity
                                               min_chunkSize = 50)

saveRDS(rcc_sim_bulk,file = 'sim_bulk.RDS')
