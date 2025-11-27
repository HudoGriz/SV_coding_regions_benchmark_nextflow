/*
========================================================================================
    SIMULATE AND BENCHMARK WORKFLOW
========================================================================================
    Simulate target regions similar to exome+UTR and benchmark SV calls against them
----------------------------------------------------------------------------------------
*/

include { SIMULATE_TARGETS } from '../modules/local/simulate_targets'
include { TRUVARI_BENCH } from '../modules/nf-core/truvari/bench/main'

workflow SIMULATE_AND_BENCHMARK {
    take:
    ch_fasta                // channel: reference FASTA
    ch_fasta_fai            // channel: reference FAI index
    ch_benchmark_vcf        // channel: truth VCF
    ch_benchmark_vcf_tbi    // channel: truth VCF index
    ch_vcfs                 // channel: [meta, vcf, tbi] - all SV caller VCFs to benchmark
    num_simulations         // val: number of simulations to run
    ch_wes_utr_targets      // channel: WES+UTR targets BED file
    ch_high_confidence_targets  // channel: high confidence regions BED

    main:
    //
    // Simulate target regions
    //
    SIMULATE_TARGETS(
        num_simulations,
        ch_fasta,
        ch_wes_utr_targets,
        ch_high_confidence_targets
    )

    //
    // Flatten simulated BED files into individual channels
    //
    ch_simulated_beds = SIMULATE_TARGETS.out.simulated_beds
        .flatten()
        .map { bed_file ->
            def bed_name = bed_file.baseName
            [
                [id: bed_name, target_set: 'simulated'],
                bed_file
            ]
        }

    //
    // Filter out Illumina WES from simulation benchmarking
    // WES data is already restricted to exome regions, so benchmarking against
    // simulated exome-like regions would be redundant
    //
    ch_vcfs_for_simulation = ch_vcfs
        .filter { meta, vcf, tbi ->
            meta.technology != 'Illumina_WES'
        }

    //
    // Combine VCFs with simulated target BEDs for benchmarking
    //
    ch_bench_input = ch_vcfs_for_simulation
        .combine(ch_simulated_beds)
        .map { vcf_meta, vcf, tbi, bed_meta, bed ->
            def combined_meta = vcf_meta + bed_meta
            combined_meta.id = "${vcf_meta.id}_${bed_meta.id}"
            [combined_meta, vcf, tbi, bed]
        }

    //
    // Prepare input for TRUVARI_BENCH
    // The module expects: tuple val(meta), path(vcf), path(tbi), path(truth_vcf), path(truth_tbi), path(bed)
    //                     tuple val(meta2), path(fasta)
    //                     tuple val(meta3), path(fai)
    //
    
    // Combine benchmark VCF with FASTA to create inputs for all samples
    ch_truvari_input = ch_bench_input
        .combine(ch_benchmark_vcf)
        .combine(ch_benchmark_vcf_tbi)
        .map { meta, vcf, tbi, bed, truth_vcf, truth_tbi ->
            [meta, vcf, tbi, truth_vcf, truth_tbi, bed]
        }

    // Prepare FASTA channel (without meta since module expects just path)
    ch_fasta_input = ch_fasta
        .map { fasta -> [[id: 'reference'], fasta] }
    
    // Prepare FAI channel (without meta since module expects just path)
    ch_fai_input = ch_fasta_fai
        .map { fai -> [[id: 'reference_index'], fai] }

    //
    // Run Truvari benchmarking on simulated targets
    //
    TRUVARI_BENCH(
        ch_truvari_input,
        ch_fasta_input,
        ch_fai_input
    )

    emit:
    simulated_beds = SIMULATE_TARGETS.out.simulated_beds
    truvari_results = TRUVARI_BENCH.out.summary
}
