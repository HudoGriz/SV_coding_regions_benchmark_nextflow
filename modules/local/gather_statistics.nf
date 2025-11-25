process GATHER_STATISTICS {
    tag "gathering statistics and generating plots"
    label 'process_low'
    publishDir "${params.outdir}/statistics", mode: 'copy'

    input:
    path truvari_results
    path run_dir

    output:
    path "plots/*", emit: plots
    path "tables/*", emit: tables
    path "summary_statistics.txt", emit: summary

    when:
    params.gather_statistics

    script:
    """
    mkdir -p plots tables
    
    # Run R script for generating plots and statistics
    Rscript ${projectDir}/bin/R/paper_plots.R \\
        ${run_dir} \\
        plots \\
        tables
    
    # Generate summary statistics
    Rscript ${projectDir}/bin/R/general_statistics.R \\
        ${run_dir} \\
        summary_statistics.txt
    """
}
