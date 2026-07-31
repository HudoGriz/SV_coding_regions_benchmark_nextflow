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
    
    // Breakend handling is a scored variable, not a fixed setting: the GIAB
    // truth sets contain no BND records, so retaining them charges every
    // breakend a caller emits as a false positive. Each target is therefore
    // benchmarked once per mode in params.bnd_modes.
    ch_bnd_modes = Channel.fromList(
        params.bnd_modes.toString().split(',').collect { it.trim() }.findAll { it }
    )

    // Create combinations of VCFs, target sets and breakend modes
    ch_benchmark_input = ch_vcfs
        .combine(ch_targets)
        .combine(ch_bnd_modes)
        .combine(ch_benchmark_vcf)
        .combine(ch_benchmark_vcf_tbi)
        .map { meta, vcf, vcf_tbi, target_name, target_bed, bnd_mode, truth_vcf, truth_tbi ->
            // Determine if WES-specific parameters should be used
            def is_wes = meta.technology == 'Illumina_WES'
            def refdist = is_wes ? params.truvari_wes_refdist : params.truvari_refdist
            def pctsize = is_wes ? params.truvari_wes_pctsize : params.truvari_pctsize
            def pctovl = is_wes ? params.truvari_wes_pctovl : params.truvari_pctovl
            def pctseq = is_wes ? params.truvari_wes_pctseq : params.truvari_pctseq

            // Add target, breakend mode and truvari parameters to metadata
            def meta_with_args = meta + [
                target: target_name,
                bnd_mode: bnd_mode,
                truvari_args: "--refdist ${refdist} --pctsize ${pctsize} --pctovl ${pctovl} --pctseq ${pctseq}"
            ]
            meta_with_args.id = "${meta.id}_${target_name}_${bnd_mode}"

            [meta_with_args, vcf, vcf_tbi, truth_vcf, truth_tbi, target_bed]
        }
    
    // Prepare reference FASTA channel for TRUVARI_BENCH
    ch_fasta_with_meta = ch_fasta.map { fasta -> [[id: 'reference'], fasta] }
    ch_fasta_fai_with_meta = ch_fasta_fai.map { fai -> [[id: 'reference'], fai] }
    
    // Run Truvari benchmarking
    TRUVARI_BENCH(
        ch_benchmark_input,
        ch_fasta_with_meta,
        ch_fasta_fai_with_meta
    )
    
    emit:
    summary = TRUVARI_BENCH.out.summary  // channel: [meta, summary.json]
}
