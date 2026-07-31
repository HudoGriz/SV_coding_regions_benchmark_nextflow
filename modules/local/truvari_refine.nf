process TRUVARI_REFINE {
    tag "$meta.id"
    label 'process_medium'

    container "${ params.truvari_container }"

    input:
    tuple val(meta), path(tp_base), path(tp_base_tbi), path(tp_comp), path(tp_comp_tbi),
          path(fn), path(fn_tbi), path(fp), path(fp_tbi), path(summary), path(params_json)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.refine.variant_summary.json"), emit: variant_summary
    tuple val(meta), path("*.refine.region_summary.json") , emit: region_summary
    tuple val(meta), path("*.refine.regions.txt")         , emit: regions
    tuple val(meta), path("*.refine.log.txt")             , emit: log
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # `truvari refine` operates on a bench output directory. TRUVARI_BENCH
    # flattens that directory when it publishes, so rebuild it here from the
    # staged files under the names refine expects. Symlinks are enough - refine
    # reads these and writes its own output alongside them.
    mkdir -p bench
    ln -s ../${tp_base}     bench/tp-base.vcf.gz
    ln -s ../${tp_base_tbi} bench/tp-base.vcf.gz.tbi
    ln -s ../${tp_comp}     bench/tp-comp.vcf.gz
    ln -s ../${tp_comp_tbi} bench/tp-comp.vcf.gz.tbi
    ln -s ../${fn}          bench/fn.vcf.gz
    ln -s ../${fn_tbi}      bench/fn.vcf.gz.tbi
    ln -s ../${fp}          bench/fp.vcf.gz
    ln -s ../${fp_tbi}      bench/fp.vcf.gz.tbi
    ln -s ../${summary}     bench/summary.json
    ln -s ../${params_json} bench/params.json

    # --reference is passed explicitly rather than relying on the path recorded
    # in params.json, which points at the original run's staging directory.
    truvari refine \\
        --reference ${fasta} \\
        --threads ${task.cpus} \\
        ${args} \\
        bench

    mv bench/refine.variant_summary.json ./${prefix}.refine.variant_summary.json
    mv bench/refine.region_summary.json  ./${prefix}.refine.region_summary.json
    mv bench/refine.regions.txt          ./${prefix}.refine.regions.txt
    mv bench/refine.log.txt              ./${prefix}.refine.log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        truvari: \$(echo \$(truvari version 2>&1) | sed 's/^Truvari v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '{}' > ${prefix}.refine.variant_summary.json
    echo '{}' > ${prefix}.refine.region_summary.json
    touch ${prefix}.refine.regions.txt
    touch ${prefix}.refine.log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        truvari: \$(echo \$(truvari version 2>&1) | sed 's/^Truvari v//')
    END_VERSIONS
    """
}
