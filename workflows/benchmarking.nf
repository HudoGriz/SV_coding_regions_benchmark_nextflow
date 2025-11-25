/*
========================================================================================
    BENCHMARKING WORKFLOW
========================================================================================
    Benchmarks SV calls against truth set using Truvari
----------------------------------------------------------------------------------------
*/

include { TRUVARI_BENCH } from '../modules/nf-core/truvari/bench/main'

workflow BENCHMARKING {
    take:
    ch_vcfs                 // channel: [meta, vcf, tbi] - VCF files to benchmark
    ch_benchmark_vcf        // channel: truth VCF file
    ch_benchmark_vcf_tbi    // channel: truth VCF index
    ch_targets              // channel: [target_name, bed] - target regions
    ch_fasta                // channel: reference FASTA
    ch_fasta_fai            // channel: reference FAI index
    
    main:
    
    // Create combinations of VCFs and target sets for benchmarking
    ch_benchmark_input = ch_vcfs
        .combine(ch_targets)
        .map { meta, vcf, vcf_tbi, target_name, target_bed ->
            // Determine if WES-specific parameters should be used
            def is_wes = meta.technology == 'Illumina_WES'
            def refdist = is_wes ? params.truvari_wes_refdist : params.truvari_refdist
            def pctsize = is_wes ? params.truvari_wes_pctsize : params.truvari_pctsize
            def pctovl = is_wes ? params.truvari_wes_pctovl : params.truvari_pctovl
            def pctseq = is_wes ? params.truvari_wes_pctseq : params.truvari_pctseq
            
            // Add target and truvari parameters to metadata
            def meta_with_args = meta + [
                target: target_name,
                truvari_args: "--refdist ${refdist} --pctsize ${pctsize} --pctovl ${pctovl} --pctseq ${pctseq}"
            ]
            
            [meta_with_args, vcf, vcf_tbi, target_bed]
        }
    
    // Prepare truth VCF channel for TRUVARI_BENCH
    ch_truth_vcf = ch_benchmark_vcf
        .combine(ch_benchmark_vcf_tbi)
        .map { vcf, tbi -> [[id: 'truth'], vcf, tbi] }
    
    // Prepare reference FASTA channel for TRUVARI_BENCH
    ch_fasta_with_index = ch_fasta
        .combine(ch_fasta_fai)
        .map { fasta, fai -> [[id: 'reference'], fasta, fai] }
    
    // Run Truvari benchmarking
    TRUVARI_BENCH(
        ch_benchmark_input,
        ch_truth_vcf,
        ch_fasta_with_index
    )
    
    emit:
    summary = TRUVARI_BENCH.out.summary  // channel: [meta, summary.json]
}
