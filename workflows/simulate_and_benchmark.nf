/*
========================================================================================
    SIMULATE AND BENCHMARK WORKFLOW
========================================================================================
    Simulate target regions similar to exome+UTR and benchmark SV calls against them
----------------------------------------------------------------------------------------
*/

include { SIMULATE_TARGETS } from '../modules/local/simulate_targets'
include { TRUVARI_BENCH } from '../modules/local/truvari'

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
    // Run Truvari benchmarking on simulated targets
    //
    TRUVARI_BENCH(
        ch_bench_input,
        ch_benchmark_vcf,
        ch_benchmark_vcf_tbi,
        ch_fasta,
        ch_fasta_fai
    )

    emit:
    simulated_beds = SIMULATE_TARGETS.out.simulated_beds
    truvari_results = TRUVARI_BENCH.out.summary
}
