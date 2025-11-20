process PBSV {
    tag "${meta.id}_${meta.tool}"
    label 'process_high'
    
    publishDir "${params.outdir}/calls/${meta.technology}/sv/pbsv", mode: 'copy'
    
    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_fai
    
    output:
    tuple val(meta), path("HG002_hs37d5.vcf.gz"), path("HG002_hs37d5.vcf.gz.tbi"), emit: vcf
    path "HG002.svsig.gz", emit: svsig
    
    script:
    """
    # Step 1: Discover structural variant signatures
    pbsv discover \\
        ${bam} \\
        HG002.svsig.gz
    
    # Step 2: Call structural variants
    pbsv call \\
        ${fasta} \\
        HG002.svsig.gz \\
        HG002_hs37d5.vcf \\
        -j ${task.cpus}
    
    # Compress and index VCF
    bgzip HG002_hs37d5.vcf
    tabix HG002_hs37d5.vcf.gz
    """
}
