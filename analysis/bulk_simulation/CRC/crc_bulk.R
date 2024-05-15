
source('../../../scripts/funcs.R')
source('../../../scripts/libraries.R')
library(Seurat)

crc_sim_bulk = list()

##################### tcag #####################
CRC_refPhi = readRDS('../../../refPhi/CRC_refPhi.RDS')
tcga_coad = deconvBenchmarking::build_tcga_obj('COAD',to_use_genes = rownames(CRC_refPhi@phi.cs))
tcga_read = deconvBenchmarking::build_tcga_obj('READ',to_use_genes = rownames(CRC_refPhi@phi.cs))

crc_sim_bulk$TCGA$simulated_bulk = cbind(tcga_coad$simulated_bulk,tcga_read$simulated_bulk)
crc_sim_bulk$TCGA$simulated_frac = rbind(tcga_coad$simulated_frac,tcga_read$simulated_frac)


#################### heter/pseudobulk from Data_Lee2020_Colorectal ##############
dataset_path = '../../../../../curated_dataset/CRC/Data_Lee2020_Colorectal/'

scExpr = Seurat::ReadMtx(mtx = paste0(dataset_path,'Exp_data_UMIcounts.mtx'),
                         cells = paste0(dataset_path,'Cells.csv'),
                         features = paste0(dataset_path,'Genes.txt'),
                         feature.column = 1,skip.cell = 1,cell.sep = ',') 
scMeta = read.delim(paste0(dataset_path,'Cells.csv'),sep =',')
scMeta = scMeta %>% column_to_rownames('cell_name')


# remove cells with no annotations/ also remove 'Mast' since there's only 1 mast cell
keep_id = which(scMeta$cell_type!='' & scMeta$cell_type!='Mast')
scExpr = scExpr[,keep_id]
scMeta = scMeta[keep_id,]

n_pseudo = length(unique(scMeta$sample))
n_simu = 50-n_pseudo

simulated_frac = fracSimulator_Dirichlet(table(scMeta$cell_type),n=n_simu,dispersion_par = 0.0005,min.frac = 0.01)
crc_sim_bulk$Lee2020 = create_pseudobulk_obj(scExpr,scMeta,unit = 'UMI',add_heter = T,simulated_frac = simulated_frac)

saveRDS(crc_sim_bulk,file = 'sim_bulk.RDS')
