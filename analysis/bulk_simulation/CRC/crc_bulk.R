
source('../../../scripts/funcs.R')
source('../../../scripts/libraries.R')
library(Seurat)

crc_sim_bulk = list()

##################### tcag #####################
CRC_refPhi = readRDS('../../../refPhi/CRC_refPhi.RDS')
crc_sim_bulk$TCGA = deconvBenchmarking::build_tcga_obj('CRC',
                                                        to_use_genes = rownames(CRC_refPhi@phi.cs))

#################### heter/pseudobulk from Data_Lee2020_Colorectal ##############
scExpr = Seurat::ReadMtx(mtx = '../../../curated_dataset/CRC/Data_Lee2020_Colorectal/Exp_data_UMIcounts.mtx',
                         cells = '../../../curated_dataset/CRC/Data_Lee2020_Colorectal/Cells.csv',
                         features = '../../../curated_dataset/CRC/Data_Lee2020_Colorectal/Genes.txt',
                         feature.column = 1,skip.cell = 1,cell.sep = ',') # less memory usage compared with readMM()
scMeta = read.delim('../../../curated_dataset/CRC/Data_Lee2020_Colorectal/Cells.csv',sep =',')
scMeta = scMeta %>% column_to_rownames('cell_name')

ct_table = table(scMeta$cell_type)
ct_table = ct_table[names(ct_table)!='']

simulated_frac = fracSimulator_Dirichlet(ct_table,n=30,dispersion_par = 0.0005,min.frac = 0.01)
crc_sim_bulk$Lee2020 = create_pseudobulk_obj(scExpr,scMeta,unit = 'UMI',add_heter = T,simulated_frac = simulated_frac)


saveRDS(crc_sim_bulk,file = 'sim_bulk.RDS')
