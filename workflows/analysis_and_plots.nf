/*
========================================================================================
    ANALYSIS AND PLOTS WORKFLOW
========================================================================================
    Gather statistics and generate plots from benchmarking results
----------------------------------------------------------------------------------------
*/

include { GATHER_STATISTICS } from '../modules/local/gather_statistics'

workflow ANALYSIS_AND_PLOTS {
    take:
    ch_truvari_results      // channel: Truvari benchmark results [meta, summary.json]

    main:
    //
    // Organize Truvari results for staging
    // We need to stage the actual result directories from work directories, not rely on publishing
    //
    ch_organized_results = ch_truvari_results
        .map { meta, summary ->
            // Get the parent directory containing all Truvari outputs
            def result_dir = summary.parent
            // Determine if this is real intervals or simulated
            def target_set = meta.target_set ?: 'real'
            // Create staging path based on target_set
            def stage_path = target_set == 'simulated' ? 
                "simulated_intervals/benchmarks/${meta.technology}/${meta.tool}/${meta.target}" :
                "real_intervals/${meta.technology}/truvari/${meta.tool}/${meta.target}"
            
            // Return a map with metadata for easier processing
            [stage_path: stage_path, result_dir: result_dir]
        }
        .collect()

    //
    // Generate statistics and plots
    //
    GATHER_STATISTICS(
        ch_organized_results
    )

    emit:
    plots = GATHER_STATISTICS.out.plots
    tables = GATHER_STATISTICS.out.tables
    // summary = GATHER_STATISTICS.out.summary  // Commented out - not currently generated
}
