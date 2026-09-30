// Drops calls genotyped homozygous reference (0/0) before any benchmark: the
// caller is stating that the sample does not carry them, and the truth sets
// count only ALT-carrying records. Calls without a genotype (./.) are kept.
// The script runs from the repository inside the analysis container, like the
// other Python steps. `removed` is the number of records dropped, so the caller
// of this module can keep using the original VCF when nothing was removed.

process EXCLUDE_HOMREF_CALLS {
    tag "${meta.technology}:${meta.tool}"
    label 'process_single'

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${prefix}.benchmarked.vcf.gz"), path("${prefix}.benchmarked.vcf.gz.tbi"), env(removed), emit: vcf
    path "${prefix}.homref_excluded.tsv", emit: counts

    script:
    prefix = task.ext.prefix ?: "${meta.technology}-${meta.tool}"
    """
    python3 ${projectDir}/bin/python/exclude_homref_calls.py \\
        --input ${vcf} \\
        --output ${prefix}.benchmarked.vcf.gz \\
        --counts ${prefix}.homref_excluded.tsv \\
        --removed-count removed_count.txt \\
        --technology ${meta.technology} \\
        --caller ${meta.tool}
    removed=\$(cat removed_count.txt)
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.technology}-${meta.tool}"
    """
    touch ${prefix}.benchmarked.vcf.gz ${prefix}.benchmarked.vcf.gz.tbi
    printf 'technology\\tcaller\\trecords\\thomref_removed\\thomref_removed_pass\\tkept\\n' > ${prefix}.homref_excluded.tsv
    removed=0
    """
}
