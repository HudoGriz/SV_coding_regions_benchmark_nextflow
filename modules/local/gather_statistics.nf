process GATHER_STATISTICS {
    tag "gathering statistics and generating plots"
    label 'process_low'

    input:
    val results_info  // List of maps with stage_path and result_dir

    output:
    path "plots/*", emit: plots, optional: true
    path "tables/*", emit: tables, optional: true

    when:
    params.gather_statistics

    script:
    // Create staging instructions for symlinks
    def staging_instructions = results_info.collect { info ->
        "mkdir -p \"\$(dirname '${info.stage_path}')\" && ln -sf '${info.result_dir}' '${info.stage_path}'"
    }.join(' && ')
    
    """
    #!/bin/bash
    set -e
    
    echo "=== Starting Analysis and Plots Generation ==="
    echo "Processing ${results_info.size()} Truvari benchmark results"
    
    # Create directory structure that R scripts expect
    mkdir -p real_intervals simulated_intervals/benchmarks plots tables
    
    # Stage all Truvari result directories using symlinks
    ${staging_instructions}
    
    echo "Directory structure created and data staged"
    
    # Run R script for generating plots and statistics
    Rscript ${projectDir}/bin/R/paper_plots.R || {
        echo "Warning: paper_plots.R encountered an error"
        echo "Check if real_intervals or simulated_intervals data exists"
    }
    
    echo "Analysis complete"
    """
}
