# Table of contents:
1)	Introduction
2)	Scripts
3)	Input and output data preparation and organization
# 1) Introduction:
The purpose of the code in this repository is to use the InstaPrism R package (https://github.com/humengying0907/InstaPrism/tree/master) to run Bayes Prism deconvolution on the following bulk RNA seq and microarray datasets of high grade serous ovarian carcinoma (HGSOC) samples:
- “SchildkrautB” – bulk RNA sequencing of HGSOC from Black patients.
- “SchildkrautW” – bulk RNA sequencing of HGSOC from White patients.
- “TCGA_bulk” – bulk RNA sequencing of HGSOC, accessed via the curatedOvarianData R package (dataset name in package: “TCGA.RNASeqV2_eset”)
- “TCGA_microarray” – microarray of HGSOC, accessed via the curatedOvarianData R package (dataset name in package: “TCGA_eset”)
- “Tothill” – microarray of HGSOC, accessed via the curatedOvarianData R package (dataset name in package: “GSE9891_eset”)
- “Yoshihara” – microarray of HGSOC, accessed via the curatedOvarianData R package (dataset name in package: “GSE32062.GPL6480_eset”)

Bayes Prism requires a single cell reference dataset to perform deconvolution. However, it has been previously shown that certain cell types, notably adipocytes, are present in bulk tumor samples but largely absent from single cell RNA sequencing results (https://www.biorxiv.org/content/10.1101/2024.04.25.590992v1). However, adipocytes can be captured using single nucleus RNA sequencing. Here, we incorporate single nucleus RNA sequencing data from adipocytes, in addition to single cell RNA sequencing data of HGSOC, in deconvolution of bulk HGSOC RNA sequencing and microarray data. 
# 2) Scripts:
Prior to running these scripts, please download the required data as outlined below in “Input and output data preparation and organization.” Please also ensure that the renv folder has been downloaded and is in the same directory as the scripts. These scripts are intended to be run in the following order:
## 1_get_data.R
This script reads the bulk RNA sequencing and microarray datasets and filters them to only include genes present in one common gene mapping list. It transforms the microarray data using 2^(...) to match the scale of the bulk RNA sequencing data values  – this is used for InstaPrism deconvolution. It also transforms the bulk RNA sequencing data using log10(...+1) to match the scale of the microarray data – this is used for clustering. All of these matrices are saved in a uniform format containing on sample ID (rows) and genes (columns); file names are appended with either “asImported” or “transformed.” It also saves any metadata information (ex. prior clustering) about the samples in a separate file for reference.
## 2_get_clustering.R
This script performs k-means clustering (k=2,3,4), NMF clustering (k=2,3,4), and consensusOV subtyping (https://github.com/bhklab/consensusOV) for each bulk dataset. For k-means clustering, it uses log10(...+1) transformed data (for RNAseq), and raw log2 data (for microarray). For NMF clustering and consensusOV subtyping, it uses raw counts (for RNAseq), and 2^(...) scaled "pseudocounts" data (for microarray). It saves the results in one csv file per dataset, where each row corresponds to one sample. It also saves a csv file containing the results for all datasets.
## 3_prepare_reference_data.R
This script creates a single-cell/single-nucleus reference matrix for use as input into InstaPrism. It requires single cell RNA sequencing data of HGSOC (from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217517) and cell type labels for those cells (from https://github.com/greenelab/deconvolution_pilot/tree/main/data/cell_labels). It also requires single nucleus RNA sequencing data of adipocytes (from http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176171). Please refer to the “Input and output data preparation and organization” subsection of README for specific files necessary. It reads in the HGSOC single cell and adipocyte single nucleus RNA sequencing data and performs some pre-processing on the adipocyte data (removing duplicate samples, removing non-adipocyte samples, removing samples with many mitochondrial gene reads, and using Seurat to remove low-quality nuclei, empty droplets, and nuclei doublets/multiplets). It then combines the single cell and single nucleus data to generate an expression matrix where each row corresponds to a gene (GeneCards symbols) and each column corresponds to a sample (each assigned a unique numerical ID). It also generates a cell type file which serves as a key, and associates each sample ID to its cell type. 
## 4_run_instaprism.R
This script runs InstaPrism, an R package that performs deconvolution with a similar but faster method to BayesPrism. First, it loads in the reference single cell plus single nucleus RNA sequencing reference dataset created in the previous script. Since the reference dataset is so large, it randomly selects 500 of each cell type to use. It then generates and saves two reference objects for input into InstaPrism, one with adipocytes and one without adipocytes. It runs InstaPrism twice on each of the six bulk datasets, both with and without the adipocytes in the reference data. Of note, Instaprism requires non-log-transformed bulk data, so it is performed on the original, non-transformed data for the bulk RNA sequencing datasets and on the 2^(…) transformed data for the microarray datasets.
## 5_visualize_instaprism_outputs.R
This script creates several figures to visualize the deconvolution results generated in the previous script, and compare the results when run with and without adipocyte single nucleus RNA sequencing data in the reference data. It creates 100% stacked bar charts to visualize the total cell proportions per bulk dataset, 100% stacked bar charts showing cell proportions per sample in each dataset, and bar charts showing the absolute change of cell type proportions in total per dataset.
# 3) Input and output data preparation and organization:
Prior to running these scripts, please ensure that the below raw data files have been downloaded and are present in a folder entitled “input_data” inside the same project and directory as the scripts. (The results will be created and stored in another folder inside the same directory entitled “output_data”.) Sample directory contents:

![Screenshot 2025-02-25 at 2 28 25 PM](https://github.com/user-attachments/assets/b44cfc2a-8242-4f5e-ba2b-f75ff6712de9)
## SchildkrautB, SchildkrautW, and reference gene list:
### From https://github.com/greenelab/hgsc_characterization:
- /reference_data/ensembl_hgnc_entrez.tsv
- /reference_data/main_AA_metadata_table.tsv
- /reference_data/main_white_metadata_table.tsv
- /data/way_pipeline_results_10removed_NeoRemoved_inclWhites/1.DataInclusion-Data-Genes/GlobalMAD_genelist.csv
- /data/way_pipeline_results_10removed_NeoRemoved_inclWhites/2.Clustering_DiffExprs-Tables-ClusterMembership/FullClusterMembership.csv
### From https://github.com/nrosed/hgsc_characterization/tree/master: 
- /data/rna_seq_pilot_and_new/salmon_normalized_filtered_for_way_pipeline.tsv
- /data/rna_seq_whites/salmon_normalized_filtered_for_way_pipeline_whites.tsv
## Single cell RNA sequencing data and cell type IDs:
### From https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217517:
#### GSM6720925 Dissociated cells, single cell sequenced, rep1
- GSM6720925_single_cell_barcodes_2251.tsv.gz
- GSM6720925_single_cell_features_2251.tsv.gz
- GSM6720925_single_cell_matrix_2251.mtx.gz
#### GSM6720926 Dissociated cells, single cell sequenced, rep2
- GSM6720926_single_cell_barcodes_2267.tsv.gz
- GSM6720926_single_cell_features_2267.tsv.gz
- GSM6720926_single_cell_matrix_2267.mtx.gz
#### GSM6720927 Dissociated cells, single cell sequenced, rep3
- GSM6720927_single_cell_barcodes_2283.tsv.gz
- GSM6720927_single_cell_features_2283.tsv.gz
- GSM6720927_single_cell_matrix_2283.mtx.gz
#### GSM6720928 Dissociated cells, single cell sequenced, rep4
- GSM6720928_single_cell_barcodes_2293.tsv.gz
- GSM6720928_single_cell_features_2293.tsv.gz
- GSM6720928_single_cell_matrix_2293.mtx.gz
#### GSM6720929 Dissociated cells, single cell sequenced, rep5
- GSM6720929_single_cell_barcodes_2380.tsv.gz
- GSM6720929_single_cell_features_2380.tsv.gz
- GSM6720929_single_cell_matrix_2380.mtx.gz
#### GSM6720930 Dissociated cells, single cell sequenced, rep6
- GSM6720930_single_cell_barcodes_2428.tsv.gz
- GSM6720930_single_cell_features_2428.tsv.gz
- GSM6720930_single_cell_matrix_2428.mtx.gz
#### GSM6720931 Dissociated cells, single cell sequenced, rep7
- GSM6720931_single_cell_barcodes_2467.tsv.gz
- GSM6720931_single_cell_features_2467.tsv.gz
- GSM6720931_single_cell_matrix_2467.mtx.gz
#### GSM6720932 Dissociated cells, single cell sequenced, rep8
- GSM6720932_single_cell_barcodes_2497.tsv.gz
- GSM6720932_single_cell_features_2497.tsv.gz
- GSM6720932_single_cell_matrix_2497.mtx.gz
### From https://github.com/greenelab/deconvolution_pilot/tree/main/data/cell_labels:
- 2251_labels.txt
- 2267_labels.txt
- 2283_labels.txt
- 2293_labels.txt
- 2380_labels.txt
- 2428_labels.txt
- 2467_labels.txt
- 2497_labels.txt
## Single nucleus adipose reference data:
### From http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176171:
#### GSM5359325 Human visceral white adipose 01-1
- GSM5359325_Hs_OAT_01-1.dge.tsv.gz
#### GSM5359326 Human visceral white adipose 01-2
- GSM5359326_Hs_OAT_01-2.dge.tsv.gz
#### GSM5359327 Human visceral white adipose 253-1
- GSM5359327_Hs_OAT_253-1.dge.tsv.gz
#### GSM5359328 Human visceral white adipose 254-1
- GSM5359328_Hs_OAT_254-1.dge.tsv.gz
#### GSM5359329 Human visceral white adipose 255-1
- GSM5359329_Hs_OAT_255-1.dge.tsv.gz
#### GSM5359330 Human visceral white adipose 256-1
- GSM5359330_Hs_OAT_256-1.dge.tsv.gz
#### GSM5359331 Human subcutaneous white adipose 01-1
- GSM5359331_Hs_SAT_01-1.dge.tsv.gz
#### GSM5359332 Human subcutaneous white adipose 02-1
- GSM5359332_Hs_SAT_02-1.dge.tsv.gz
#### GSM5359333 Human subcutaneous white adipose 04-1
- GSM5359333_Hs_SAT_04-1.dge.tsv.gz
#### GSM5359334 Human subcutaneous white adipose 253-1
- GSM5359334_Hs_SAT_253-1.dge.tsv.gz
#### GSM5359335 Human subcutaneous white adipose 254-1
- GSM5359335_Hs_SAT_254-1.dge.tsv.gz
#### GSM5359336 Human subcutaneous white adipose 255-1
- GSM5359336_Hs_SAT_255-1.dge.tsv.gz
#### GSM5359337 Human subcutaneous white adipose 256-1
- GSM5359337_Hs_SAT_256-1.dge.tsv.gz
#### GSM5820679	Human visceral white adipose 09-1
- GSM5820679_Hs_OAT_09-1.dge.tsv.gz
#### GSM5820680	Human visceral white adipose 10-1
- GSM5820680_Hs_OAT_10-1.dge.tsv.gz
#### GSM5820684	Human subcutaneous white adipose 09-1
- GSM5820684_Hs_SAT_09-1.dge.tsv.gz
#### GSM5820685	Human subcutaneous white adipose 10-1
- GSM5820685_Hs_SAT_10-1.dge.tsv.gz
#### GSM5820686	Human subcutaneous white adipose 11-1
- GSM5820686_Hs_SAT_11-1.dge.tsv.gz
#### Metadata
- GSE176171_cell_metadata.tsv.gz
