process DOWNLOAD_BAM {
    tag "${meta.id}"
    label 'process_low'
    
    publishDir "${output_dir}", mode: 'copy'
    
    input:
    tuple val(meta), val(url), val(output_dir)
    
    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam
    
    script:
    def bam_filename = url.tokenize('/').last()
    """
    # Download BAM file
    wget -O ${bam_filename} ${url}
    
    # Download BAI index (should be at same location with .bai extension)
    wget -O ${bam_filename}.bai ${url}.bai
    """
    
    stub:
    def bam_filename = url.tokenize('/').last()
    """
    touch ${bam_filename}
    touch ${bam_filename}.bai
    """
}

// Alias processes for different technologies/genome builds
process DOWNLOAD_ILLUMINA_WES_BAM {
    tag "${meta.id}"
    publishDir "${output_dir}", mode: 'copy'
    
    input:
    tuple val(meta), val(url), val(output_dir)
    
    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam
    
    script:
    def bam_filename = url.tokenize('/').last()
    """
    wget -O ${bam_filename} ${url}
    wget -O ${bam_filename}.bai ${url}.bai
    """
}

process DOWNLOAD_ILLUMINA_WGS_BAM {
    tag "${meta.id}"
    publishDir "${output_dir}", mode: 'copy'
    
    input:
    tuple val(meta), val(url), val(output_dir)
    
    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam
    
    script:
    def bam_filename = url.tokenize('/').last()
    """
    wget -O ${bam_filename}.bai ${url}.bai
    wget -O ${bam_filename} ${url}
    """
}

process DOWNLOAD_ILLUMINA_WGS_BAM_GRCH38 {
    tag "${meta.id}"
    publishDir "${output_dir}", mode: 'copy'
    
    input:
    tuple val(meta), val(url), val(output_dir)
    
    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam
    
    script:
    def bam_filename = url.tokenize('/').last()
    """
    wget -O ${bam_filename}.bai ${url}.bai
    wget -O ${bam_filename} ${url}
    """
}

process DOWNLOAD_PACBIO_BAM {
    tag "${meta.id}"
    publishDir "${output_dir}", mode: 'copy'
    
    input:
    tuple val(meta), val(url), val(output_dir)
    
    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam
    
    script:
    def bam_filename = url.tokenize('/').last()
    """
    wget -O ${bam_filename} ${url}
    wget -O ${bam_filename}.bai ${url}.bai
    """
}

process DOWNLOAD_PACBIO_BAM_GRCH38 {
    tag "${meta.id}"
    publishDir "${output_dir}", mode: 'copy'
    
    input:
    tuple val(meta), val(url), val(output_dir)
    
    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam
    
    script:
    def bam_filename = url.tokenize('/').last()
    """
    wget -O ${bam_filename}.bai ${url}.bai
    wget -O ${bam_filename} ${url}
    """
}

process DOWNLOAD_ONT_BAM {
    tag "${meta.id}"
    publishDir "${output_dir}", mode: 'copy'
    
    input:
    tuple val(meta), val(url), val(output_dir)
    
    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam
    
    script:
    def bam_filename = url.tokenize('/').last()
    """
    wget -O ${bam_filename} ${url}
    wget -O ${bam_filename}.bai ${url}.bai
    """
}

process DOWNLOAD_ONT_BAM_GRCH38 {
    tag "${meta.id}"
    publishDir "${output_dir}", mode: 'copy'
    
    input:
    tuple val(meta), val(url), val(output_dir)
    
    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam
    
    script:
    def bam_filename = url.tokenize('/').last()
    """
    wget -O ${bam_filename}.bai ${url}.bai
    wget -O ${bam_filename} ${url}
    """
}
