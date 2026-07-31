process TRUVARI_REFINE {
    tag "$meta.id"
    label 'process_medium'

    container "${ params.truvari_container }"

    input:
    tuple val(meta), path(tp_base), path(tp_base_tbi), path(tp_comp), path(tp_comp_tbi),
          path(fn), path(fn_tbi), path(fp), path(fp_tbi), path(summary), path(params_json),
          path(candidates)
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
    # refine's default --regions is <benchdir>/candidate.refine.bed
    ln -s ../${candidates}  bench/candidate.refine.bed

    # --reference is passed explicitly rather than relying on the path in
    # params.json, which points at the original run's staging directory.
    #
    # When no region holds both an unmatched call and something to harmonize it
    # against, refine logs "No regions to be refined", writes no summary, and
    # still exits 0. Illumina WES on the gene panel does exactly that: 4 TP,
    # 0 FP, 1,278 FN, so nothing can be rescued. Branch on whether the output
    # actually exists rather than on the exit code, which does not distinguish
    # the two cases.
    set +e
    truvari refine \\
        --reference ${fasta} \\
        --threads ${task.cpus} \\
        ${args} \\
        bench > refine.console.txt 2>&1
    rc=\$?
    set -e
    cat refine.console.txt

    if [ -f bench/refine.variant_summary.json ]; then
        mv bench/refine.variant_summary.json ./${prefix}.refine.variant_summary.json
        mv bench/refine.region_summary.json  ./${prefix}.refine.region_summary.json
        mv bench/refine.regions.txt          ./${prefix}.refine.regions.txt
        mv bench/refine.log.txt              ./${prefix}.refine.log.txt
    elif grep -q "No regions to be refined" refine.console.txt; then
        echo "${prefix}: nothing to harmonize; carrying the unrefined result forward"
        cp ${summary} ./${prefix}.refine.variant_summary.json
        echo '{"note": "no regions to be refined", "refined_regions": 0}' \\
            > ./${prefix}.refine.region_summary.json
        printf 'chrom\\tstart\\tend\\trefined\\n' > ./${prefix}.refine.regions.txt
        cp refine.console.txt ./${prefix}.refine.log.txt
    else
        echo "${prefix}: truvari refine produced no summary (exit \$rc)" >&2
        # never exit 0 here: reaching this branch means refine neither wrote a
        # result nor reported that there was nothing to do
        [ \$rc -ne 0 ] && exit \$rc
        exit 1
    fi

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
