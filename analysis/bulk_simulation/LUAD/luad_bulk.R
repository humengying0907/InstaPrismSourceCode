
source('../../../scripts/funcs.R')
source('../../../scripts/libraries.R')
library(Seurat)

luad_sim_bulk = list()

############################ TCGA ############################
LUAD_refPhi = readRDS('../../../refPhi/LUAD_refPhi.RDS')
luad_sim_bulk$TCGA = deconvBenchmarking:::build_tcga_obj('LUAD',
                                                         to_use_genes = rownames(LUAD_refPhi@phi.cs))

#################### pseudobulk plus heter from Data_Kim2020_Lung ############
dataset_path = '../../../../../curated_dataset/LUAD/Data_Kim2020_Lung/'


scExpr = Seurat::ReadMtx(mtx = paste0(dataset_path,'Exp_data_UMIcounts.mtx'),
                         cells = paste0(dataset_path,'Cells.csv'),
                         features = paste0(dataset_path,'Genes.txt'),
                         feature.column = 1,skip.cell = 1,cell.sep = ',') 

scMeta = read.delim(paste0(dataset_path,'Cells.csv'),sep = ',') %>% column_to_rownames('cell_name')

ct_table = table(scMeta$cell_type)
ct_table = ct_table[names(ct_table)!='']

n_pseudo = length(unique(scMeta$sample))
n_simu = 50-n_pseudo

simulated_frac = fracSimulator_Dirichlet(ct_table,n=n_simu,dispersion_par = 0.0005,min.frac = 0.01)

luad_sim_bulk$Kim2020 = create_pseudobulk_obj(scExpr,scMeta,unit = 'UMI',min.cells = 3,
                                              add_heter = T,simulated_frac = simulated_frac)

################# pseudobulk from Qian2020 ################## 
dataset_path = '../../../../../curated_dataset/LUAD/Data_Qian2020_Lung/'

scExpr = Seurat::ReadMtx(mtx = paste0(dataset_path,'Exp_data_UMIcounts.mtx'),
                         cells = paste0(dataset_path,'Cells.csv'),
                         features = paste0(dataset_path,'Genes.txt'),
                         feature.column = 1,skip.cell = 1,cell.sep = ',')
scMeta = read.delim(paste0(dataset_path,'Cells.csv'),sep = ',') %>% column_to_rownames('cell_name')

# only use LUAD samples
to_keep = which(scMeta$disease == 'LUAD' & scMeta$cell_type!='')


scExpr = scExpr[,to_keep]
scMeta = scMeta[to_keep,]

ct_table = table(scMeta$cell_type)
n_pseudo = length(unique(scMeta$sample))
n_simu = 50-n_pseudo

simulated_frac = fracSimulator_Dirichlet(ct_table,n=n_simu,dispersion_par = 0.0005,min.frac = 0.01)

luad_sim_bulk$Qian2020 = create_pseudobulk_obj(scExpr = scExpr,
                                               scMeta = scMeta,
                                               colnames_of_cellType = 'cell_type',
                                               colnames_of_sample = 'sample',
                                               unit = 'UMI',
                                               add_heter = T,
                                               simulated_frac = simulated_frac,
                                               min_chunkSize = 40)

saveRDS(luad_sim_bulk,file = 'sim_bulk.RDS')  
