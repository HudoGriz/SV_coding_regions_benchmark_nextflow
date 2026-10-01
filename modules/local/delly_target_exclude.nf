/*
    Delly exclude list for capture data: the release exclude template plus every
    base outside the capture targets. Delly's authors describe it as designed for
    WGS, with limited applicability to exome data, and advise keeping calls whose
    breakpoints lie in targeted sequence (dellytools/delly issues 207 and 275).
    Excluding everything else makes Delly's discovery region the capture targets,
    the same regions Manta uses as call regions on the WES data.
*/
process DELLY_TARGET_EXCLUDE {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(targets)
    path fai
    path template

    output:
    tuple val(meta), path("${prefix}.delly_exclude.tsv"), emit: exclude
    tuple val(meta), path("${prefix}.delly_exclude.counts.tsv"), emit: counts

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def read_targets = targets.name.endsWith('.gz') ? "gzip -dc ${targets}" : "cat ${targets}"
    def read_template = template ? "cat ${template}" : ":"
    """
    ${read_targets} | awk -F'\\t' 'BEGIN { OFS = "\\t" } !/^(#|track|browser)/ && NF >= 3 { print \$1, \$2, \$3 }' \\
        | LC_ALL=C sort -k1,1 -k2,2n > targets.sorted.bed

    # Complement of the targets over every contig of the reference.
    awk -F'\\t' 'BEGIN { OFS = "\\t" }
        function close_contig() { if (pos < len[chrom]) print chrom, pos, len[chrom]; done[chrom] = 1 }
        FNR == NR { len[\$1] = \$2; order[++n] = \$1; next }
        !(\$1 in len) { next }
        \$1 != chrom { if (chrom != "") close_contig(); chrom = \$1; pos = 0 }
        { if (\$2 > pos) print chrom, pos, \$2; if (\$3 > pos) pos = \$3 }
        END {
            if (chrom != "") close_contig()
            for (i = 1; i <= n; i++) if (!(order[i] in done)) print order[i], 0, len[order[i]]
        }' ${fai} targets.sorted.bed > outside_targets.bed

    { ${read_template}; cat outside_targets.bed; } > ${prefix}.delly_exclude.tsv

    awk -F'\\t' 'BEGIN { OFS = "\\t"; print "set", "intervals", "bp" }
        FNR == NR { tb += \$3 - \$2; tn++; next }
        { ob += \$3 - \$2; on++ }
        END { print "targets", tn, tb; print "outside_targets", on, ob }' targets.sorted.bed outside_targets.bed \\
        > ${prefix}.delly_exclude.counts.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.delly_exclude.tsv ${prefix}.delly_exclude.counts.tsv
    """
}
