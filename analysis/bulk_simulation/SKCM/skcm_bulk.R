
source('../../../scripts/funcs.R')
source('../../../scripts/libraries.R')
library(Seurat)

SKCM_refPhi = readRDS('../../../../refPhi_cs/SKCM_refPhi.RDS')

skcm_sim_bulk = list()

####################### TCGA ######################
skcm_sim_bulk$TCGA = deconvBenchmarking:::build_tcga_obj('SKCM',
                                                         to_use_genes = rownames(SKCM_refPhi@phi.cs),
                                                         purity_methods =  c('ESTIMATE', 'ABSOLUTE', 'LUMP', 'IHC', 'CPE','ABSOLUTE_GDC'))

################### pseudo plus heter from Data_Jerby-Arnon2018_Skin ###############
dataset_path = '../../../../../curated_dataset/SKCM/Data_Jerby-Arnon2018_Skin/'

scExpr = Seurat::ReadMtx(mtx = paste0(dataset_path,'Exp_data_TPM.mtx'),
                         cells = paste0(dataset_path,'Cells.csv'),
                         features = paste0(dataset_path,'Genes.txt'),
                         feature.column = 1,skip.cell = 1,cell.sep = ',') 

scMeta = read.delim(paste0(dataset_path,'Cells.csv'),sep = ',') %>% column_to_rownames('cell_name')

# only use data from new cohort and cells with cell-type annotation
to_keep = rownames(scMeta)[scMeta$cell_type!='' & scMeta$cell_cohort == 'New']
scExpr = scExpr[,to_keep]
scMeta = scMeta[to_keep,]

ct_table = table(scMeta$cell_type)
ct_table = ct_table[names(ct_table)!='']
n_pseudo = length(unique(scMeta$sample))
n_simu = 50-n_pseudo


simulated_frac = fracSimulator_Dirichlet(ct_table,n=n_simu,dispersion_par = 0.0005,min.frac = 0.01)
skcm_sim_bulk$Jerby2018 = create_pseudobulk_obj(scExpr,scMeta,unit = 'UMI',min.cells = 3,
                                                add_heter = T,simulated_frac = simulated_frac)

saveRDS(skcm_sim_bulk,file = 'sim_bulk.RDS')
