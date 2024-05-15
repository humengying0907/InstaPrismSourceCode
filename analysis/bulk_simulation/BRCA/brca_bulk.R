
source('../../../scripts/funcs.R')
source('../../../scripts/libraries.R')
library(Seurat)

brca_sim_bulk = list()

##################### tcag #####################
BRCA_refPhi = readRDS('../../../refPhi/BRCA_refPhi.RDS')
brca_sim_bulk$TCGA = deconvBenchmarking::build_tcga_obj('BRCA',
                                                         to_use_genes = rownames(BRCA_refPhi@phi.cs))


############ pseudobulk plus heter simulated bulk samples from Qian2020 ##########
dataset_path = '../../../../../curated_dataset/BRCA/Data_Qian2020_Breast/'

scExpr = Seurat::ReadMtx(mtx = paste0(dataset_path,'Exp_data_UMIcounts.mtx'),
                         cells = paste0(dataset_path,'Cells.csv'),
                         features = paste0(dataset_path,'Genes.txt'),
                         feature.column = 1,skip.cell = 1,cell.sep = ',') 
scMeta = read.delim(paste0(dataset_path,'Cells.csv'),sep =',')
scMeta = scMeta %>% column_to_rownames('cell_name')

n_pseudo = length(unique(scMeta$sample))
n_simu = 50-n_pseudo

# remove cells with no annotations
keep_id = which(scMeta$cell_type!='')
scExpr = scExpr[,keep_id]
scMeta = scMeta[keep_id,]

simulated_frac = fracSimulator_Dirichlet(table(scMeta$cell_type),n=n_simu,dispersion_par = 0.0005,min.frac = 0.01)
brca_sim_bulk$Qian2020 = create_pseudobulk_obj(scExpr,scMeta,unit = 'UMI',add_heter = T,
                                               simulated_frac = simulated_frac)


########## pseudobulk plus heter simulated bulk samples from BRCA_Bassez2021_immunotherapy ##########
sce = readRDS('../../../../../curated_dataset/BRCA/BRCA_Bassez2021_immunotherapy/cohort1/sce.RDS')
scMeta = as.data.frame(colData(sce))

n_pseudo = length(unique(scMeta$patient_id))
n_simu = 50-n_pseudo
simulated_frac = fracSimulator_Dirichlet(table(scMeta$cell_type),n=n_simu,dispersion_par = 0.0005,min.frac = 0.01)
brca_sim_bulk$Bassez2021 = create_pseudobulk_obj(sce@assays@data@listData$cpm,
                                                 scMeta,
                                                 colnames_of_sample = 'patient_id',
                                                 colnames_of_cellType = 'cell_type',
                                                 unit = 'UMI',add_heter = T,
                                                 simulated_frac = simulated_frac)

saveRDS(brca_sim_bulk,file = 'sim_bulk.RDS')



