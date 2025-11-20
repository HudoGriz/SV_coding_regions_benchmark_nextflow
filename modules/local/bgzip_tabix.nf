process BGZIP_TABIX {
    tag "${meta.id}"
    label 'process_low'
    
    container 'community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa'
    
    input:
    tuple val(meta), path(vcf)
    
    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.tool}"
    """
    # Compress VCF with bgzip
    bgzip -c ${vcf} > ${prefix}.vcf.gz
    
    # Index with tabix
    tabix -p vcf ${prefix}.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bgzip: \$(bgzip --version 2>&1 | head -1 | sed 's/^.*bgzip //' | sed 's/ .*\$//')
        tabix: \$(tabix --version 2>&1 | head -1 | sed 's/^.*tabix //' | sed 's/ .*\$//')
    END_VERSIONS
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.tool}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bgzip: \$(bgzip --version 2>&1 | head -1 | sed 's/^.*bgzip //' | sed 's/ .*\$//')
        tabix: \$(tabix --version 2>&1 | head -1 | sed 's/^.*tabix //' | sed 's/ .*\$//')
    END_VERSIONS
    """
}
