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
    // Extract just the summary.json files from the [meta, summary.json] tuples
    // then collect all results into a single list
    //
    ch_collected_results = ch_truvari_results
        .map { meta, summary -> summary }
        .collect()

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
