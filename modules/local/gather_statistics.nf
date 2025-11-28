process GATHER_STATISTICS {
    tag "gathering statistics and generating plots"
    label 'process_low'

    input:
    val ready  // Signal that truvari results are ready
    path run_dir

    output:
    path "plots/*", emit: plots, optional: true
    path "tables/*", emit: tables, optional: true
    path "summary_statistics.txt", emit: summary

    when:
    params.gather_statistics

    script:
    """
    #!/bin/bash
    set -e
    
    echo "=== Starting Analysis and Plots Generation ==="
    echo "Run directory: ${run_dir}"
    echo "Working directory: \$(pwd)"
    echo "Processing ${ready} Truvari benchmark results"
    echo ""
    
    # Create output directories
    mkdir -p plots tables
    
    echo "=== Running paper_plots.R ==="
    # Run R script for generating plots and statistics
    Rscript ${projectDir}/bin/R/paper_plots.R \\
        ${run_dir} \\
        plots \\
        tables || {
        echo "Warning: paper_plots.R encountered an error, but continuing..."
        echo "Check if real_intervals or simulated_intervals data exists in ${run_dir}"
    }
    
    echo ""
    echo "=== Running general_statistics.R ==="
    # Generate summary statistics
    Rscript ${projectDir}/bin/R/general_statistics.R \\
        ${run_dir} \\
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
