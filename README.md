# InstaPrismSourceCode

## Built-in Reference for InstaPrism

| reference name | tumor type                       | #cells used for reference construction | #cell types/cell states | umap | citation | download |
|----------------|----------------------------------|----------------------------------------|-------------------------|------|----------|----------|
| BRCA_refPhi    | breast cancer                    | 100,064                                    | 8/76                             |[UMAP](https://singlecell.broadinstitute.org/single_cell/study/SCP1039/a-single-cell-and-spatially-resolved-atlas-of-human-breast-cancers)      | [Wu et al. 2021](https://www.nature.com/articles/s41588-021-00911-1) | [↓](https://github.com/humengying0907/InstaPrismSourceCode/raw/main/refPhi/BRCA_refPhi.RDS) |
| CRC_refPhi     | colorectal cancer                | 371,223                                    | 15/98                           |[UMAP](https://singlecell.broadinstitute.org/single_cell/study/SCP1162/human-colon-cancer-atlas-c295?scpbr=human-cell-atlas-main-collection)      | [Pelka et al. 2021](https://www.cell.com/cell/fulltext/S0092-8674(21)00945-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421009454%3Fshowall%3Dtrue) |  [↓](https://github.com/humengying0907/InstaPrismSourceCode/raw/main/refPhi/CRC_refPhi.RDS) |
| GBM_refPhi     | glioblastoma                     | 338,564                                    | 10/57                             |[cellxgeneLink](https://cellxgene.cziscience.com/collections/999f2a15-3d7e-440b-96ae-2c806799c08c), [UMAP](https://cellxgene.cziscience.com/e/56c4912d-2bae-4b64-98f2-af8a84389208.cxg/)      | [Ruiz et al. 2022](https://www.biorxiv.org/content/10.1101/2022.08.27.505439v1) |[↓](https://github.com/humengying0907/InstaPrismSourceCode/raw/main/refPhi/GBM_refPhi.RDS) |
| LUAD_refPhi    | lung adenocarcinomas             | 118,293                                    | 13/77                            |[UMAP](https://www.weizmann.ac.il/sites/3CA/study-data/umap/20115)      | [Xing et al. 2021](https://www.science.org/doi/10.1126/sciadv.abd9738) | [↓](https://github.com/humengying0907/InstaPrismSourceCode/raw/main/refPhi/LUAD_refPhi.RDS) |
| OV_refPhi      | ovarian cancer                   | 929,690                                    | 9/40                             |[cellxgeneLink](https://cellxgene.cziscience.com/collections/4796c91c-9d8f-4692-be43-347b1727f9d8), [UMAP](https://cellxgene.cziscience.com/e/b252b015-b488-4d5c-b16e-968c13e48a2c.cxg/)      | [Vazquez et al. 2022](https://www.nature.com/articles/s41586-022-05496-1) | [↓](https://github.com/humengying0907/InstaPrismSourceCode/raw/main/refPhi/OV_refPhi.RDS) |
| RCC_refPhi     | clear cell renal cell carcinoma  | 270,855                                    | 11/106                            |[cellxgeneLink](https://cellxgene.cziscience.com/collections/f7cecffa-00b4-4560-a29a-8ad626b8ee08), [UMAP](https://cellxgene.cziscience.com/e/5af90777-6760-4003-9dba-8f945fec6fdf.cxg/)      | [Li et al. 2022](https://www.cell.com/cancer-cell/fulltext/S1535-6108(22)00548-7) | [↓](https://github.com/humengying0907/InstaPrismSourceCode/raw/main/refPhi/RCC_refPhi.RDS) |
| SKCM_refPhi    | skin cutaneous melanoma          | 4,645                                      | 8/23                             |[UMAP](https://www.weizmann.ac.il/sites/3CA/study-data/umap/20111)      | [Tirosh et al. 2016](https://www.science.org/doi/10.1126/science.aad0501) | [↓](https://github.com/humengying0907/InstaPrismSourceCode/raw/main/refPhi/SKCM_refPhi.RDS) |

## Reference validation pipeline
1. To reproduce the reference validation results from the paper or to test the reference performance with your own data, clone the repository and go into the InstaPrismSourceCode directory.
```
git clone https://github.com/humengying0907/InstaPrismSourceCode.git && cd InstaPrismSourceCode
```

2. Make sure that all the packages in `scripts/libraries.R` are installed
```
# To install InstaPrism & deconvBenchmarking package
library("devtools");
install_github("humengying0907/InstaPrism")
install_github("humengying0907/deconvBenchmarking")
```

3. Make simulated bulk samples for reference validation. Code to reproduce bulk simulation for each cancer type is listed in `./analysis/bulk_simulation/"cancer type"`. A `sim_bulk.RDS` file will be generated and used in the subsequent validation step. Single cell data used for bulk simulation is downloaded from [3CA](https://www.weizmann.ac.il/sites/3CA/) repository.

4. Reference validation. The script used for deconvolution and performance evaluation is `./analysis/refPhi_validation/evalu_pipeline.R`

```
cd analysis/refPhi_validation
Rscript evalu_pipeline.R -h
```
    ## usage: evalu_pipeline.R [-h] [--tumorType <character>]
    ##                         [--refName <character> [<character> ...]]
    ##                         [--output <character>] [--niter <integer>]
    ##                         [--ncore <integer>] [--updateReference <logical>]
    ##                         [--key <character>] [--saveDeconvRes <logical>]
    ## 
    ## InstaPrism reference evaluation pipeline
    ## 
    ## optional arguments:
    ##   -h, --help            show this help message and exit
    ##   --tumorType <character>, -t <character>
    ##                         Folder name of the bulk dataset for evaluation. The
    ##                         'sim_bulk.RDS' file located within
    ##                         '../bulk_simulation/<tumorType>' will be utilized for
    ##                         deconvolution.
    ##   --refName <character> [<character> ...], -n <character> [<character> ...]
    ##                         Name of the reference to test. The reference stored in
    ##                         '../../refPhi/<refName>_refPhi.RDS' will be loaded as
    ##                         the reference. To test multiple references, separate
    ##                         each name by a blank space.
    ##   --output <character>, -o <character>
    ##                         Output folder name. A 'performance/<output>/'
    ##                         directory will be created to store the deconvolution
    ##                         results and plots. If not specified, a
    ##                         'performance/<tumorType>/' directory will be created
    ##                         to store the results
    ##   --niter <integer>     Number of iterations for InstaPrism() function.
    ##                         [default: 400]
    ##   --ncore <integer>     Number of threads. [default: 16]
    ##   --updateReference <logical>
    ##                         A logical variable to determine whether to include
    ##                         updated reference in the evaulation. [default: FALSE]
    ##   --key <character>     Name of the malignant cell type in the reference,
    ##                         required only when updateReference = TRUE. When
    ##                         evaluating multiple references, these references need
    ##                         to have the same key.
    ##   --saveDeconvRes <logical>
    ##                         A logical variable to determine whether to save the
    ##                         deconvolution results. [default: FALSE]

Run the following commands to generate the performance summary object `theta_performance.RDS` and the associated pdf files that visualize the performance of the deconvolution results. An example output is listed in `./analysis/refPhi_validation/performance/CRC`. A summarized performance plot can be found at `./analysis/refPhi_validation/performance/performance_summary.png`

```
Rscript evalu_pipeline.R -t BRCA -n BRCA --updateReference TRUE --key "Cancer Epithelial"
Rscript evalu_pipeline.R -t CRC -n CRC --updateReference TRUE --key EpiT
Rscript evalu_pipeline.R -t GBM -n GBM --updateReference TRUE --key Malignant
Rscript evalu_pipeline.R -t LUAD -n LUAD --updateReference TRUE --key Malignant
Rscript evalu_pipeline.R -t OV -n OV --updateReference TRUE --key malignant
Rscript evalu_pipeline.R -t RCC -n RCC --updateReference TRUE --key Malignant
Rscript evalu_pipeline.R -t SKCM -n SKCM --updateReference TRUE --key "Melanoma cells"
```

## Citation
M. Hu and M. Chikina, “InstaPrism: an R package for fast implementation of BayesPrism.” bioRxiv, p. 2023.03.07.531579, Mar. 12, 2023.
doi: https://doi.org/10.1101/2023.03.07.531579
