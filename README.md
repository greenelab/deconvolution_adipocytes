# deconvolution_adipocytes
The following repository is a work in progress. The aim is to deconvolve HGSOC samples with added adipocyte signal.


## Prior to running the first script (/scripts/1_get_data_and_clustering.R) in the repository, please make sure you have the following files/data:

### From https://github.com/greenelab/hgsc_characterization:

### Reference Data:

/reference_data/ensembl_hgnc_entrez.tsv

/reference_data/main_AA_metadata_table.tsv

/reference_data/main_white_metadata_table.tsv

/reference_data/consensusOV_res_05.tsv

/reference_data/protype_predictions.tsv

/reference_data/protype_predictions_merged.tsv

### Data Files:

/data/way_pipeline_results_10removed_NeoRemoved_inclWhites/1.DataInclusion-Data-Genes/CommonGenes_genelist.csv

/data/way_pipeline_results_10removed_NeoRemoved_inclWhites/1.DataInclusion-Data-Genes/GlobalMAD_genelist.csv

/data/way_pipeline_results_10removed_NeoRemoved_inclWhites/2.Clustering_DiffExprs-Tables-ClusterMembership/FullClusterMembership.csv

### Clustering and Expression Data:

/data/mayo/MayoEset.Rda

TCGA dataset (loaded via data(TCGA_eset), part of the curatedOvarianData package)

Tothill dataset (loaded via data(GSE9891_eset), part of the curatedOvarianData package)

Yoshihara dataset (loaded via data(GSE32062.GPL6480_eset), part of the curatedOvarianData package)

### Utility Scripts (Sourced Files):

/analyze_rnaseq/plot_utils.R

/comparison/utils/file_processing_utils.R


### From https://github.com/nrosed/hgsc_characterization/tree/master:

### Data Files:

/data/rna_seq_pilot_and_new/salmon_normalized_filtered_for_way_pipeline.tsv

/data/rna_seq_pilot_and_new/salmon_normalized_filtered_for_way_pipeline_MAD.tsv

/data/rna_seq_whites/salmon_normalized_filtered_for_way_pipeline_whites.tsv

/data/rna_seq_whites/salmon_normalized_filtered_for_way_pipeline_whites_MAD.tsv


## Prior to running the third script (/scripts/3_prepare_reference_data.R) in the repository, please make sure you have the following file/data:

### From https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217517:

### GSM6720925	Dissociated cells, single cell sequenced, rep1 

GSM6720925_single_cell_barcodes_2251.tsv.gz

GSM6720925_single_cell_features_2251.tsv.gz

GSM6720925_single_cell_matrix_2251.mtx.gz

### GSM6720926	Dissociated cells, single cell sequenced, rep2

GSM6720926_single_cell_barcodes_2267.tsv.gz

GSM6720926_single_cell_features_2267.tsv.gz

GSM6720926_single_cell_matrix_2267.mtx.gz

### GSM6720927	Dissociated cells, single cell sequenced, rep3

GSM6720927_single_cell_barcodes_2283.tsv.gz

GSM6720927_single_cell_features_2283.tsv.gz

GSM6720927_single_cell_matrix_2283.mtx.gz

### GSM6720928	Dissociated cells, single cell sequenced, rep4

GSM6720928_single_cell_barcodes_2293.tsv.gz

GSM6720928_single_cell_features_2293.tsv.gz

GSM6720928_single_cell_matrix_2293.mtx.gz

### GSM6720929	Dissociated cells, single cell sequenced, rep5

GSM6720929_single_cell_barcodes_2380.tsv.gz

GSM6720929_single_cell_features_2380.tsv.gz

GSM6720929_single_cell_matrix_2380.mtx.gz

### GSM6720930	Dissociated cells, single cell sequenced, rep6

GSM6720930_single_cell_barcodes_2428.tsv.gz

GSM6720930_single_cell_features_2428.tsv.gz

GSM6720930_single_cell_matrix_2428.mtx.gz

### GSM6720931	Dissociated cells, single cell sequenced, rep7

GSM6720931_single_cell_barcodes_2467.tsv.gz

GSM6720931_single_cell_features_2467.tsv.gz

GSM6720931_single_cell_matrix_2467.mtx.gz

### GSM6720932	Dissociated cells, single cell sequenced, rep8

GSM6720932_single_cell_barcodes_2497.tsv.gz

GSM6720932_single_cell_features_2497.tsv.gz

GSM6720932_single_cell_matrix_2497.mtx.gz

### From https://github.com/greenelab/deconvolution_pilot/tree/main/data/cell_labels:

2251_labels.txt

2267_labels.txt

2283_labels.txt

2293_labels.txt

2380_labels.txt

2428_labels.txt

2467_labels.txt

2497_labels.txt





