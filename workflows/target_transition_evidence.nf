/*
========================================================================================
    TARGET TRANSITION EVIDENCE WORKFLOW
========================================================================================
*/

include { TARGET_TRANSITION_AUDIT; MERGE_TARGET_TRANSITION_EVIDENCE; PLOT_TARGET_TRANSITION_EVIDENCE; SIMULATION_TRANSITION_AUDIT; MERGE_SIMULATION_TRANSITION_EVIDENCE } from '../modules/local/target_transition_evidence'

workflow TARGET_TRANSITION_EVIDENCE {
    take:
    ch_benchmarks       // channel: [meta, tp-base, tp-comp, fn, fp]
    ch_targets          // channel: [target_name, bed]
    ch_simulation_benchmarks // channel: [meta, tp-base, tp-comp, fn, fp]
    ch_simulated_beds   // channel: collection of simulation BED files

    main:
    ch_hci = ch_benchmarks
        .filter { meta, tp_base, tp_comp, fn, fp ->
            meta.target == params.transition_hci_target &&
                (params.transition_include_wes || meta.technology != 'Illumina_WES')
        }
        .map { meta, tp_base, tp_comp, fn, fp ->
            def key = "${meta.technology}:${meta.tool}:${meta.bnd_mode ?: 'default'}"
            [key, meta, [tp_base, tp_comp, fn, fp]]
        }

    ch_restricted = ch_benchmarks
        .filter { meta, tp_base, tp_comp, fn, fp ->
            meta.target == params.transition_target &&
                (params.transition_include_wes || meta.technology != 'Illumina_WES')
        }
        .map { meta, tp_base, tp_comp, fn, fp ->
            def key = "${meta.technology}:${meta.tool}:${meta.bnd_mode ?: 'default'}"
            [key, meta, [tp_base, tp_comp, fn, fp]]
        }

    ch_target_bed = ch_targets
        .filter { target_name, bed -> target_name == params.transition_target }
        .map { target_name, bed -> bed }

    ch_pairs = ch_hci
        .join(ch_restricted)
        .combine(ch_target_bed)
        .map { key, hci_meta, hci_files, target_meta, target_files, target_bed ->
            [key, hci_meta, hci_files, target_meta, target_files, target_bed]
        }

    TARGET_TRANSITION_AUDIT(ch_pairs)

    ch_all_transitions = TARGET_TRANSITION_AUDIT.out.transitions.collect()
    MERGE_TARGET_TRANSITION_EVIDENCE(ch_all_transitions)

    ch_merged_transitions = MERGE_TARGET_TRANSITION_EVIDENCE.out.evidence
        .flatten()
        .filter { file -> file.name.endsWith('.transitions.tsv') }
    PLOT_TARGET_TRANSITION_EVIDENCE(ch_merged_transitions)

    // Materialize the BED lookup as a value channel. A keyed queue-to-queue
    // join/combine consumes the single BED emission for a simulation and can
    // therefore pair it with only one caller. The value is broadcast to every
    // caller/simulation pair, and the matching BED is selected from the map.
    ch_sim_beds_by_id = ch_simulated_beds
        .flatten()
        .map { bed -> [bed.baseName, bed] }
        .collect(flat: false)
        .map { pairs -> pairs.collectEntries { simulation_id, bed -> [(simulation_id): bed] } }

    ch_simulation_inputs = ch_simulation_benchmarks
        .map { meta, tp_base, tp_comp, fn, fp ->
            def pipeline_key = "${meta.technology}:${meta.tool}:${meta.bnd_mode ?: 'default'}"
            [pipeline_key, meta.target, meta, [tp_base, tp_comp, fn, fp]]
        }

    ch_simulation_pairs = ch_hci
        .combine(ch_simulation_inputs, by: 0)
        .map { pipeline_key, hci_meta, hci_files, simulation_id, simulation_meta, simulation_files ->
            [simulation_id, "${pipeline_key}:${simulation_id}", hci_meta, hci_files, simulation_meta, simulation_files]
        }
        .combine(ch_sim_beds_by_id)
        .map { simulation_id, pair_id, hci_meta, hci_files, simulation_meta, simulation_files, bed_by_id ->
            def bed = bed_by_id[simulation_id]
            if (bed == null) {
                error("No simulated BED found for ${simulation_id}")
            }
            [pair_id, hci_meta, hci_files, simulation_meta, simulation_files, bed]
        }

    SIMULATION_TRANSITION_AUDIT(ch_simulation_pairs)
    ch_all_simulation_transitions = SIMULATION_TRANSITION_AUDIT.out.transitions.collect()
    MERGE_SIMULATION_TRANSITION_EVIDENCE(ch_all_simulation_transitions)

    emit:
    evidence = MERGE_TARGET_TRANSITION_EVIDENCE.out.evidence
    figures = PLOT_TARGET_TRANSITION_EVIDENCE.out.figures
    simulation_evidence = MERGE_SIMULATION_TRANSITION_EVIDENCE.out.evidence
}
