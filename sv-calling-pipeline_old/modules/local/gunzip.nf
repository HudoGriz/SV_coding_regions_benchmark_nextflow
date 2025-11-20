process GUNZIP {
    tag "${archive.name}"
    label 'process_low'
    
    publishDir "${publishDir}", mode: 'copy', enabled: params.publish_gunzip
    
    input:
    path archive
    val publishDir
    
    output:
    path "${output_name}", emit: file
    
    script:
    output_name = archive.name.replaceAll(/\.gz$/, '')
    """
    gunzip -c ${archive} > ${output_name}
    """
    
    stub:
    output_name = archive.name.replaceAll(/\.gz$/, '')
    """
    touch ${output_name}
    """
}
