process DOWNLOAD_TANDEM_REPEATS {
    tag "${meta.id}"
    label 'process_low'
    
    publishDir "${output_dir}", mode: 'copy'
    
    input:
    tuple val(meta), val(url), val(output_dir)
    
    output:
    path "*.bed", emit: bed
    
    script:
    """
    wget -O human_hs37d5.trf.bed ${url}
    """
    
    stub:
    """
    touch human_hs37d5.trf.bed
    """
}

process DOWNLOAD_TANDEM_REPEATS_GRCH38 {
    tag "${meta.id}"
    label 'process_low'
    
    publishDir "${output_dir}", mode: 'copy'
    
    input:
    tuple val(meta), val(url), val(output_dir)
    
    output:
    path "*.bed", emit: bed
    
    script:
    """
    wget -O human_GRCh38_no_alt_analysis_set.trf.bed ${url}
    """
    
    stub:
    """
    touch human_GRCh38_no_alt_analysis_set.trf.bed
    """
}

process DOWNLOAD_GENCODE_GTF {
    tag "${meta.id}"
    label 'process_low'
    
    publishDir "${output_dir}", mode: 'copy'
    
    input:
    tuple val(meta), val(url), val(output_dir)
    
    output:
    path "*.gtf.gz", emit: gtf
    
    script:
    """
    wget -O gencode.v19.annotation.gtf.gz ${url}
    """
    
    stub:
    """
    touch gencode.v19.annotation.gtf.gz
    """
}

process DOWNLOAD_GENCODE_GTF_GRCH38 {
    tag "${meta.id}"
    label 'process_low'
    
    publishDir "${output_dir}", mode: 'copy'
    
    input:
    tuple val(meta), val(url), val(output_dir)
    
    output:
    path "*.gtf.gz", emit: gtf
    
    script:
    """
    wget -O gencode.v49.annotation.gtf.gz ${url}
    """
    
    stub:
    """
    touch gencode.v49.annotation.gtf.gz
    """
}
