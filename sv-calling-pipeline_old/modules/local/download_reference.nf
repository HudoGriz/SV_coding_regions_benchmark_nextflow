process DOWNLOAD_REFERENCE_GENOME {
    tag "${meta.id}"
    label 'process_low'
    
    publishDir "${output_dir}", mode: 'copy'
    
    input:
    tuple val(meta), val(url), val(output_dir)
    
    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    
    script:
    """
    # Download and decompress reference genome
    wget -O human_hs37d5.fasta.gz ${url}
    gunzip human_hs37d5.fasta.gz
    """
    
    stub:
    """
    touch human_hs37d5.fasta
    """
}

process DOWNLOAD_REFERENCE_GENOME_GRCH38 {
    tag "${meta.id}"
    label 'process_low'
    
    publishDir "${output_dir}", mode: 'copy'
    
    input:
    tuple val(meta), val(url), val(output_dir)
    
    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    
    script:
    """
    # Download and decompress GRCh38 reference genome
    wget -O human_GRCh38_no_alt_analysis_set.fasta.gz ${url}
    gunzip human_GRCh38_no_alt_analysis_set.fasta.gz
    """
    
    stub:
    """
    touch human_GRCh38_no_alt_analysis_set.fasta
    """
}
