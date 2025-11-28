process SIMULATE_TARGETS {
    tag "simulating ${num_simulations} target regions"
    label 'process_high'
    publishDir "${params.outdir}/simulations/target_regions", mode: 'copy'
    
    input:
    val num_simulations
    path reference_fasta
    path wes_utr_targets
    path high_confidence_targets

    output:
    path "simulated_targets/*.bed", emit: simulated_beds
    path "simulated_targets/", emit: simulation_dir

    script:
    """
    echo "Creating simulated_targets directory..."
    mkdir -p simulated_targets || { echo "Failed to create directory!"; exit 1; }
    
    echo "Running R script..."
    Rscript ${projectDir}/bin/R/simulate_targets.R \\
        ${num_simulations} \\
        simulated_targets \\
        ${wes_utr_targets} \\
        ${high_confidence_targets}
    """
}
