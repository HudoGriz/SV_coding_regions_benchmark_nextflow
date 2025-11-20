process SAMTOOLS_FAIDX {
    tag "${fasta.name}"
    label 'process_low'
    
    container 'file:///singularity_images/samtools_latest.sif'
    
    publishDir "${params.references_dir}", mode: 'copy'
    
    input:
    path fasta
    
    output:
    path "${fasta}", emit: fasta
    path "${fasta}.fai", emit: fai
    
    script:
    """
    samtools faidx ${fasta}
    """
    
    stub:
    """
    touch ${fasta}.fai
    """
}
