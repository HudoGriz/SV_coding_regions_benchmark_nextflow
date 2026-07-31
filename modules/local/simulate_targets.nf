process SIMULATE_TARGETS {
    tag "simulating ${num_simulations} target regions"
    label 'process_high'
    
    input:
    val num_simulations
    path reference_fasta
    path wes_utr_targets
    path high_confidence_targets

    output:
    path "simulated_targets/*.bed", emit: simulated_beds
    path "simulated_targets/", emit: simulation_dir

    script:
    // task.cpus MUST be passed explicitly. The R script used to call
    // detectCores(), which inside the container reports the whole compute node
    // rather than the cgroup allocation; that forked 255 workers into a 64 GB
    // request and the out-of-memory handler killed most of them, silently
    // reducing a 500-set empirical null to 111 sets.
    """
    echo "Creating simulated_targets directory..."
    mkdir -p simulated_targets || { echo "Failed to create directory!"; exit 1; }

    echo "Running R script with ${task.cpus} worker(s)..."
    Rscript ${projectDir}/bin/R/simulate_targets.R \\
        ${num_simulations} \\
        simulated_targets \\
        ${wes_utr_targets} \\
        ${high_confidence_targets} \\
        ${task.cpus}

    echo "Simulated interval sets written: \$(ls -1 simulated_targets/*.bed | wc -l)"
    """
}
