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

