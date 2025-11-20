process BEDTOOLS_INTERSECT {
    tag "${meta.id}"
    label 'process_low'
    
    container 'file:///singularity_images/bedtools_latest.sif'
    
    publishDir "${params.references_dir}", mode: 'copy'
    
    input:
    tuple val(meta), path(bed_a), path(bed_b)
    
    output:
    tuple val(meta), path("*.bed"), emit: bed
    
    script:
    def output_name = meta.output_name ?: "${bed_a.baseName}.intersect.bed"
    """
    bedtools intersect -a ${bed_a} -b ${bed_b} > ${output_name}
    """
    
    stub:
    def output_name = meta.output_name ?: "${bed_a.baseName}.intersect.bed"
    """
    touch ${output_name}
    """
}
