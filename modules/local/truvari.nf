process TRUVARI_BENCH {
    tag "${meta.technology}_${meta.tool}_${meta.target}"
    label 'process_medium'
    
    publishDir "${params.outdir}/real_intervals/${meta.technology}/truvari_benchmark/${meta.tool}/${meta.target}", 
        mode: 'copy'
    
    input:
    tuple val(meta), path(call_vcf), path(call_vcf_tbi), path(target_bed)
    path benchmark_vcf
    path benchmark_vcf_tbi
    path fasta
    path fasta_fai
    
    output:
    tuple val(meta), path("truvari_output"), emit: results
    path "truvari_output/summary.json", emit: summary
    
    script:
    // Determine if WES-specific parameters should be used
    def is_wes = meta.technology == 'Illumina_WES'
    def refdist = is_wes ? params.truvari_wes_refdist : params.truvari_refdist
    def pctsize = is_wes ? params.truvari_wes_pctsize : params.truvari_pctsize
    def pctovl = is_wes ? params.truvari_wes_pctovl : params.truvari_pctovl
    def pctseq = is_wes ? params.truvari_wes_pctseq : params.truvari_pctseq
    
    """
    # Run Truvari benchmarking
    truvari bench \\
        -b ${benchmark_vcf} \\
        -c ${call_vcf} \\
        -f ${fasta} \\
        --includebed ${target_bed} \\
        -o truvari_output \\
        --dup-to-ins \\
        --passonly \\
        --refdist ${refdist} \\
        --pctsize ${pctsize} \\
        --pctovl ${pctovl} \\
        --pctseq ${pctseq}
    """
}
