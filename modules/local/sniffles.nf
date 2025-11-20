process SNIFFLES {
    tag "${meta.id}_${meta.tool}"
    label 'process_high'
    
    publishDir "${params.outdir}/calls/${meta.technology}/sv/sniffles", mode: 'copy'
    
    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_fai
    path tandem_repeats
    
    output:
    tuple val(meta), path("HG002_hs37d5.vcf.gz"), path("HG002_hs37d5.vcf.gz.tbi"), emit: vcf
    
    script:
    def tandem_repeats_arg = tandem_repeats.name != 'NO_FILE' ? "--tandem-repeats ${tandem_repeats}" : ""
    """
    # Run Sniffles
    sniffles \\
        --input ${bam} \\
        --reference ${fasta} \\
        ${tandem_repeats_arg} \\
        --threads ${task.cpus} \\
        --vcf HG002_hs37d5.vcf
    
    # Compress and index VCF
    bgzip HG002_hs37d5.vcf
    tabix HG002_hs37d5.vcf.gz
    """
}
