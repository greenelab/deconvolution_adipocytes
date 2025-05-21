#!/usr/bin/env nextflow
/*
 * main.nf — run 5 R scripts in strict order
 *
 * You already version-lock packages with renv.lock.
 * Each script therefore starts with `renv::load()` (or restore),
 * so the only thing we need is the R interpreter on PATH.
 */

workflow {

    /*
     * 1) Download / tidy raw data
     */
    process GET_DATA {
        tag '1_get_data'
        script:
        """
        cd ${params.projectDir}
        Rscript --vanilla ${params.scriptDir}/1_get_data.R
        """
    }

    /*
     * 2) Cluster single-cell data
     */
    process GET_CLUSTERING {
        tag '2_get_clustering'
        script:
        """
        cd ${params.projectDir}
        Rscript --vanilla ${params.scriptDir}/2_get_clustering.R
        """
    }
    GET_CLUSTERING.after GET_DATA        // enforce order

    /*
     * 3) Build reference matrices
     */
    process PREP_REF_DATA {
        tag '3_prepare_reference_data'
        script:
        """
        cd ${params.projectDir}
        Rscript --vanilla ${params.scriptDir}/3_prepare_reference_data.R
        """
    }
    PREP_REF_DATA.after GET_CLUSTERING

    /*
     * 4) Deconvolution with InstaPrism
     */
    process RUN_INSTAPRISM {
        tag '4_run_instaprism'
        script:
        """
        cd ${params.projectDir}
        Rscript --vanilla ${params.scriptDir}/4_run_instaprism.R
        """
    }
    RUN_INSTAPRISM.after PREP_REF_DATA

    /*
     * 5) Visualisation step
     */
    process VISUALISE {
        tag '5_visualize'
        script:
        """
        cd ${params.projectDir}
        Rscript --vanilla ${params.scriptDir}/5_visualize_instaprism_outputs.R
        """
    }
    VISUALISE.after RUN_INSTAPRISM
}
