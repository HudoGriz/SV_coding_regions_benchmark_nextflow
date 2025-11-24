process TABIX_VCF {
    tag "${vcf.name}"
    label 'process_low'
    publishDir "${params.references_dir}", mode: 'copy'
    
    input:
    path vcf
    
    output:
    path "${vcf}", emit: vcf
    path "${vcf}.tbi", emit: tbi
    
    script:
    """
    tabix -p vcf ${vcf}
    """
    
    stub:
    """
    touch ${vcf}.tbi
    """
}
