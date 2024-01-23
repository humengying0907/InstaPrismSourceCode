
source('../../../scripts/funcs.R')
source('../../../scripts/libraries.R')
library(Seurat)

luad_sim_bulk = list()

############################ TCGA ############################
LUAD_refPhi = readRDS('../../../refPhi/LUAD_refPhi.RDS')
luad_sim_bulk$TCGA = deconvBenchmarking:::build_tcga_obj('LUAD',
                                                         to_use_genes = rownames(LUAD_refPhi@phi.cs))

#################### pseudobulk plus heter from Data_Kim2020_Lung ############
# only contains 14 samples, metastatic samples not included
scExpr = Seurat::ReadMtx(mtx = '../../../../curated_dataset/LUAD/Data_Kim2020_Lung/Exp_data_UMIcounts.mtx',
                         cells = '../../../../curated_dataset/LUAD/Data_Kim2020_Lung/Cells.csv',
                         features = '../../../../curated_dataset/LUAD/Data_Kim2020_Lung/Genes.txt',
                         feature.column = 1,skip.cell = 1,cell.sep = ',') # less memory usage compared with readMM()

scMeta = read.delim('../../../../curated_dataset/LUAD/Data_Kim2020_Lung/Cells.csv',sep = ',') %>% column_to_rownames('cell_name')

ct_table = table(scMeta$cell_type)
ct_table = ct_table[names(ct_table)!='']

simulated_frac = fracSimulator_Dirichlet(ct_table,n=30,dispersion_par = 0.0005,min.frac = 0.01)

luad_sim_bulk$Kim2020 = create_pseudobulk_obj(scExpr,scMeta,unit = 'UMI',min.cells = 3,add_heter = T,simulated_frac = simulated_frac)

################# pseudobulk from Qian2020 ################## 
scExpr = Seurat::ReadMtx(mtx = '../../../../curated_dataset/LUAD/Data_Qian2020_Lung/Exp_data_UMIcounts.mtx',
                         cells = '../../../../curated_dataset/LUAD/Data_Qian2020_Lung/Cells.csv',
                         features = '../../../../curated_dataset/LUAD/Data_Qian2020_Lung/Genes.txt',
                         feature.column = 1,skip.cell = 1,cell.sep = ',') # less memory usage compared with readMM()

scMeta = read.delim('../../../../curated_dataset/LUAD/Data_Qian2020_Lung/Cells.csv',sep = ',') %>% column_to_rownames('cell_name')

# only use LUAD samples
to_keep = which(scMeta$disease == 'LUAD' & scMeta$cell_type!='')


scExpr = scExpr[,to_keep]
scMeta = scMeta[to_keep,]


ct_table = table(scMeta$cell_type)
simulated_frac = fracSimulator_Dirichlet(ct_table,n=50,dispersion_par = 0.0005,min.frac = 0.01)

luad_sim_bulk$Qian2020 = create_pseudobulk_obj(scExpr = scExpr,
                                               scMeta = scMeta,
                                               colnames_of_cellType = 'cell_type',
                                               colnames_of_sample = 'sample',
                                               unit = 'UMI',
                                               add_heter = T,
                                               simulated_frac = simulated_frac,
                                               min_chunkSize = 40)

saveRDS(luad_sim_bulk,file = 'sim_bulk.RDS')  
