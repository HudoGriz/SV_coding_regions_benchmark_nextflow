process CUTESV {
    tag "${meta.id}_${meta.tool}"
    label 'process_high'
    
    publishDir "${params.outdir}/calls/${meta.technology}/sv/cutesv", mode: 'copy'
    
    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_fai
    
    output:
    tuple val(meta), path("HG002_hs37d5.vcf.gz"), path("HG002_hs37d5.vcf.gz.tbi"), emit: vcf
    path "work_dir", emit: work_dir
    
    script:
    """
    # Create working directory for cuteSV
    mkdir -p work_dir
    
    # Run cuteSV
    cuteSV \\
        ${bam} \\
        ${fasta} \\
        HG002_hs37d5.vcf \\
        work_dir \\
        --threads ${task.cpus}
    
    # Compress and index VCF
    bgzip HG002_hs37d5.vcf
    tabix HG002_hs37d5.vcf.gz
    """
}
