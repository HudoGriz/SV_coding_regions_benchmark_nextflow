// The scripts run from the repository (bin/python) inside the analysis
// container, which supplies the pinned Python environment. That mirrors how
// the R scripts are run, and means the code that produced a result is the code
// at the commit recorded in the run manifest, not a copy baked into the image.

process TARGET_TRANSITION_AUDIT {
    tag "${target_meta.technology}:${target_meta.tool}:${target_meta.target}"
    label 'process_medium'

    input:
    tuple val(pair_id), val(hci_meta), path(hci_files), val(target_meta), path(target_files), path(target_bed)

    output:
    path "evidence/*.transitions.tsv", emit: transitions
    path "evidence/*.mechanisms.tsv", emit: mechanisms
    path "evidence/*.boundary.tsv", emit: boundaries
    path "evidence/*.metadata.json", emit: metadata

    script:
    def assembly = params.reference_assembly ?: 'unknown'
    def pipeline = "${target_meta.technology}:${target_meta.tool}"
    def safe_id = pair_id.replaceAll(/[^A-Za-z0-9_.-]/, '_')
    """
    # Record identity includes the sample genotype (target_transition_audit.py).
    mkdir -p hci target evidence
    cp ${hci_files.join(' ')} hci/
    cp ${target_files.join(' ')} target/

    python3 ${projectDir}/bin/python/target_transition_audit.py \
        --assembly '${assembly}' \
        --hci-bench hci \
        --target-bench target \
        --target-bed ${target_bed} \
        --target-name '${params.transition_target_label}' \
        --pipeline '${pipeline}' \
        --prefix evidence/${safe_id}
    """

    stub:
    def safe_id = pair_id.replaceAll(/[^A-Za-z0-9_.-]/, '_')
    """
    mkdir -p evidence
    touch evidence/${safe_id}.transitions.tsv evidence/${safe_id}.mechanisms.tsv evidence/${safe_id}.boundary.tsv
    echo '{}' > evidence/${safe_id}.metadata.json
    """
}

process MERGE_TARGET_TRANSITION_EVIDENCE {
    tag "target-transition evidence"
    label 'process_low'

    input:
    path transition_tables

    output:
    path "target_transition_evidence.*", emit: evidence

    script:
    """
    mkdir -p batches
    cp ${transition_tables.join(' ')} batches/
    python3 ${projectDir}/bin/python/merge_transition_audits.py \
        --input-dir batches \
        --prefix target_transition_evidence
    """

    stub:
    """
    touch target_transition_evidence.transitions.tsv target_transition_evidence.mechanisms.tsv
    """
}

process PLOT_TARGET_TRANSITION_EVIDENCE {
    tag "target-transition figures"
    label 'process_low'

    input:
    path merged_transitions

    output:
    path "figures/*.png", emit: figures

    script:
    """
    mkdir -p figures
    python3 ${projectDir}/bin/python/plot_target_boundary_mechanisms.py \
        --transitions ${merged_transitions} \
        --output-dir figures
    """

    stub:
    """
    mkdir -p figures
    touch figures/target_boundary_mechanisms.png
    """
}

process SIMULATION_TRANSITION_AUDIT {
    tag "${target_meta.technology}:${target_meta.tool}:${target_meta.target}"
    label 'process_low'

    input:
    tuple val(pair_id), val(hci_meta), path(hci_files), val(target_meta), path(target_files), path(target_bed)

    output:
    path "simulation_evidence/*.transitions.tsv", emit: transitions
    path "simulation_evidence/*.mechanisms.tsv", emit: mechanisms
    path "simulation_evidence/*.boundary.tsv", emit: boundaries
    path "simulation_evidence/*.metadata.json", emit: metadata

    script:
    def assembly = params.reference_assembly ?: 'unknown'
    def pipeline = "${target_meta.technology}:${target_meta.tool}"
    def safe_id = pair_id.replaceAll(/[^A-Za-z0-9_.-]/, '_')
    """
    # Record identity includes the sample genotype (target_transition_audit.py).
    mkdir -p hci target simulation_evidence
    cp ${hci_files.join(' ')} hci/
    cp ${target_files.join(' ')} target/

    python3 ${projectDir}/bin/python/target_transition_audit.py \
        --assembly '${assembly}' \
        --hci-bench hci \
        --target-bench target \
        --target-bed ${target_bed} \
        --target-name '${target_meta.target}' \
        --pipeline '${pipeline}' \
        --prefix simulation_evidence/${safe_id}
    """

    stub:
    def safe_id = pair_id.replaceAll(/[^A-Za-z0-9_.-]/, '_')
    """
    mkdir -p simulation_evidence
    touch simulation_evidence/${safe_id}.transitions.tsv simulation_evidence/${safe_id}.mechanisms.tsv simulation_evidence/${safe_id}.boundary.tsv
    echo '{}' > simulation_evidence/${safe_id}.metadata.json
    """
}

process MERGE_SIMULATION_TRANSITION_EVIDENCE {
    tag "simulation transition evidence"
    label 'process_low'

    input:
    path transition_tables

    output:
    path "simulation_transition_evidence.*", emit: evidence

    script:
    """
    mkdir -p batches
    cp ${transition_tables.join(' ')} batches/
    python3 ${projectDir}/bin/python/merge_transition_audits.py \
        --input-dir batches \
        --prefix simulation_transition_evidence
    """

    stub:
    """
    touch simulation_transition_evidence.transitions.tsv simulation_transition_evidence.mechanisms.tsv
    """
}
