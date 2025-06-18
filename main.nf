#!/usr/bin/env nextflow
/*
 * main.nf — run 6 R scripts in strict order
 *
 * You already version-lock packages with renv.lock.
 * Each script therefore starts with `renv::load()`,
 * so the only thing we need is the R interpreter on PATH.
 */

nextflow.enable.dsl = 2

/*
 * 1) Load the renv
 */
process LOAD_RENV {
    tag '01_load_renv'
    
    input:
    val dummy_input // Receive the trigger
    
    output:
    val true, emit: renv_loaded
    
    script:
    """
    cd ${params.projectDir}
    Rscript --vanilla ${params.scriptDir}/01_load_renv.R
    """
}

/*
 * 2) Download / tidy raw data
 */
process GET_DATA {
    tag '02_get_data'
    
    input:
    val dummy_input
    
    output:
    val true, emit: data_ready
    
    script:
    """
    cd ${params.projectDir}
    Rscript --vanilla ${params.scriptDir}/02_get_data.R
    """
}

/*
 * 3) Cluster single-cell data
 */
process GET_CLUSTERING {
    tag '03_get_clustering'
    
    input:
    val dummy_input
    
    output:
    val true, emit: clustering_done
    
    script:
    """
    cd ${params.projectDir}
    Rscript --vanilla ${params.scriptDir}/03_get_clustering.R
    """
}

/*
 * 4) Build reference matrices
 */
process PREP_REF_DATA {
    tag '04_prepare_reference_data'
    
    input:
    val dummy_input
    
    output:
    val true, emit: ref_data_ready
    
    script:
    """
    cd ${params.projectDir}
    Rscript --vanilla ${params.scriptDir}/04_prepare_reference_data.R
    """
}

/*
 * 5) Deconvolution with InstaPrism
 */
process RUN_INSTAPRISM {
    tag '05_run_instaprism'
    
    input:
    val dummy_input
    
    output:
    val true, emit: instaprism_complete
    
    script:
    """
    cd ${params.projectDir}
    Rscript --vanilla ${params.scriptDir}/05_run_instaprism.R
    """
}

/*
 * 6) Visualisation step
 */
process VISUALISE {
    tag '06_visualize_instaprism_outputs'
    
    input:
    val dummy_input
    
    script:
    """
    cd ${params.projectDir}
    Rscript --vanilla ${params.scriptDir}/06_visualize_instaprism_outputs.R
    """
    // No output needed as this is the last step
}

workflow {
    // Create initial trigger channel
    trigger_channel = Channel.value(true)
    
    // Run processes in strict order using the output of the previous as input for the next
    LOAD_RENV(trigger_channel)
    GET_DATA(LOAD_RENV.out.renv_loaded)
    GET_CLUSTERING(GET_DATA.out.data_ready)
    PREP_REF_DATA(GET_CLUSTERING.out.clustering_done)
    RUN_INSTAPRISM(PREP_REF_DATA.out.ref_data_ready)
    VISUALISE(RUN_INSTAPRISM.out.instaprism_complete)
}