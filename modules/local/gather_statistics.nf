process GATHER_STATISTICS {
    tag "gathering statistics and generating plots"
    label 'process_low'

    input:
    val results_info  // List of maps with stage_path and result_dir

    output:
    path "plots/*", emit: plots, optional: true
    path "tables/*", emit: tables, optional: true
    path "summary_statistics.txt", emit: summary
    path "real_intervals/**", emit: real_intervals, optional: true
    path "simulated_intervals/**", emit: simulated_intervals, optional: true

    when:
    params.gather_statistics

    script:
    // Extract paths and directories from the list of maps
    def stage_paths_list = results_info.collect { it.stage_path }
    def result_dirs_list = results_info.collect { it.result_dir }
    
    // Create a staging map for the input directive
    def staging_instructions = results_info.collect { info ->
        "mkdir -p \"\$(dirname '${info.stage_path}')\" && ln -sf '${info.result_dir}' '${info.stage_path}'"
    }.join(' && ')
    
    """
    #!/bin/bash
    set -e
    
    echo "=== Starting Analysis and Plots Generation ==="
    echo "Working directory: \$(pwd)"
    echo "Processing ${results_info.size()} Truvari benchmark results"
    echo ""
    
    # Create directory structure that R scripts expect
    echo "=== Creating directory structure ==="
    mkdir -p real_intervals simulated_intervals/benchmarks
    
    # Stage all Truvari result directories using symlinks
    echo "=== Staging Truvari result directories ==="
    ${staging_instructions}
    
    echo ""
    echo "Staging complete. Directory structure:"
    find real_intervals simulated_intervals -type l -o -type d 2>/dev/null | head -20
    echo ""
    
    # Create output directories
    mkdir -p plots tables
    
    echo "=== Running paper_plots.R ==="
    # Run R script for generating plots and statistics
    # Pass current directory as run_dir since we staged everything here
    Rscript ${projectDir}/bin/R/paper_plots.R \\
        \$(pwd) \\
        plots \\
        tables || {
        echo "Warning: paper_plots.R encountered an error, but continuing..."
        echo "Check if real_intervals or simulated_intervals data exists"
    }
    
    echo ""
    echo "=== Running general_statistics.R ==="
    # Generate summary statistics
    Rscript ${projectDir}/bin/R/general_statistics.R \\
        \$(pwd) \\
        summary_statistics.txt || {
        echo "Error: general_statistics.R failed"
        exit 1
    }
    
    echo ""
    echo "=== Analysis Complete ==="
    echo "Generated files:"
    ls -lh plots/ tables/ summary_statistics.txt 2>/dev/null || echo "Some output directories may be empty"
    """
}
