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
    ch_gencode_gtf          // channel: GENCODE GTF file
    ch_benchmark_vcf        // channel: truth VCF
    ch_benchmark_vcf_tbi    // channel: truth VCF index
    ch_vcfs                 // channel: [meta, vcf, tbi] - all SV caller VCFs to benchmark
    num_simulations         // val: number of simulations to run

    main:
    //
    // Simulate target regions
    //
    SIMULATE_TARGETS(
        num_simulations,
        ch_fasta,
        ch_gencode_gtf
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
    // Combine VCFs with simulated target BEDs for benchmarking
    //
    ch_bench_input = ch_vcfs
        .combine(ch_simulated_beds)
        .map { vcf_meta, vcf, tbi, bed_meta, bed ->
            def combined_meta = vcf_meta + bed_meta
            combined_meta.id = "${vcf_meta.id}_${bed_meta.id}"
            [combined_meta, vcf, tbi, bed]
        }

    //
    // Prepare input for TRUVARI_BENCH: combine truth VCF+index and FASTA+index into proper format
    //
    ch_truth_vcf_with_index = ch_benchmark_vcf
        .combine(ch_benchmark_vcf_tbi)
        .map { vcf, tbi -> [[id: 'truth'], vcf, tbi] }
    
    ch_fasta_with_index = ch_fasta
        .combine(ch_fasta_fai)
        .map { fasta, fai -> [[id: 'reference'], fasta, fai] }

    //
    // Remap bench input to match nf-core module expectations
    //
    ch_truvari_input = ch_bench_input
        .map { meta, vcf, tbi, bed ->
            [meta, vcf, tbi, bed]
        }

    //
    // Run Truvari benchmarking on simulated targets
    //
    TRUVARI_BENCH(
        ch_truvari_input,
        ch_truth_vcf_with_index,
        ch_fasta_with_index
    )

    emit:
    simulated_beds = SIMULATE_TARGETS.out.simulated_beds
    truvari_results = TRUVARI_BENCH.out.summary
}
