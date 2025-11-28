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
    ch_run_dir              // channel: Run directory with all results

    main:
    //
    // Wait for all Truvari results to complete, then emit a ready signal
    // We don't need to stage the files since R scripts read directly from run_dir
    //
    ch_ready = ch_truvari_results
        .map { meta, summary -> 1 }  // Convert to simple counter
        .collect()                    // Collect all results
        .map { it.size() }            // Return count as ready signal

    //
    // Generate statistics and plots
    //
    GATHER_STATISTICS(
        ch_ready,
        ch_run_dir
    )

    emit:
    plots = GATHER_STATISTICS.out.plots
    tables = GATHER_STATISTICS.out.tables
    summary = GATHER_STATISTICS.out.summary
}
