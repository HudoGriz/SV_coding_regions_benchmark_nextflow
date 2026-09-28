// Pads a target BED symmetrically, clips it to the high-confidence intervals
// and merges overlapping components: the padding procedure of the published
// boundary-sensitivity experiment. The script runs from the repository inside
// the analysis container, like the other Python steps.

process PAD_TARGET_BED {
    tag "${target_name} +${padding} bp"
    label 'process_single'

    input:
    tuple val(target_name), path(target_bed), val(padding), path(allowed_bed)

    output:
    tuple val(target_name), val(padding), path("${target_name}.pad${padding}.bed"), emit: bed

    script:
    """
    python3 ${projectDir}/bin/python/pad_target_bed.py \\
        --target ${target_bed} \\
        --allowed ${allowed_bed} \\
        --padding ${padding} \\
        --output ${target_name}.pad${padding}.bed
    """

    stub:
    """
    printf '1\\t1000\\t1200\\n' > ${target_name}.pad${padding}.bed
    """
}
