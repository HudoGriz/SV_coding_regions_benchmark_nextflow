process GATHER_STATISTICS {
    tag "gathering statistics and generating plots"
    label 'process_low'

    input:
    tuple val(stage_paths), path(result_dirs, stageAs: 'staged_*')

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
    echo "Working directory: \$(pwd)"
    echo "Processing ${stage_paths.size()} Truvari benchmark results"
    echo ""
    
    # Create directory structure that R scripts expect
    echo "=== Creating directory structure ==="
    mkdir -p real_intervals simulated_intervals/benchmarks
    
    # Stage all Truvari result directories using symlinks
    echo "=== Staging Truvari result directories ==="
    
    # Create symlinks for each result directory at the expected path
    # We need to iterate through the lists in parallel
    count=0
    for staged_dir in staged_*; do
        if [ -d "\$staged_dir" ]; then
            # Get the corresponding path from the list
            target_path=\$(echo '${stage_paths.join('\n')}' | sed -n "\$((count+1))p")
            
            # Create parent directories
            mkdir -p "\$(dirname "\$target_path")"
            
            # Create symlink
            ln -sf "\$(pwd)/\$staged_dir" "\$target_path"
            
            echo "Linked \$staged_dir -> \$target_path"
            count=\$((count+1))
        fi
    done
    
    echo ""
    echo "Directory structure:"
    find real_intervals simulated_intervals -type l -o -type d | head -20
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
