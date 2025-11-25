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
    ch_truvari_results      // channel: Truvari benchmark results
    ch_run_dir              // channel: Run directory with all results

    main:
    //
    // Collect all Truvari results
    //
    ch_collected_results = ch_truvari_results.collect()

    //
    // Generate statistics and plots
    //
    GATHER_STATISTICS(
        ch_collected_results,
        ch_run_dir
    )

    emit:
    plots = GATHER_STATISTICS.out.plots
    tables = GATHER_STATISTICS.out.tables
    summary = GATHER_STATISTICS.out.summary
}
