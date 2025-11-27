process SIMULATE_TARGETS {
    tag "simulating ${num_simulations} target regions"
    label 'process_medium'

    input:
    val num_simulations
    path reference_fasta
    path wes_utr_targets          // Add this
    path high_confidence_targets   // Add this

    output:
    path "simulated_targets/*.bed", emit: simulated_beds
    path "simulated_targets/", emit: simulation_dir

    script:
    """
    mkdir -p simulated_targets
    
    Rscript ${projectDir}/bin/R/simulate_targets.R \\
        ${num_simulations} \\
        simulated_targets \\
        ${wes_utr_targets} \\
        ${high_confidence_targets}
    """
}
