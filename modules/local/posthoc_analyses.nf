// Post-hoc analyses of the benchmarks, run as pipeline tasks when
// --posthoc_analyses is set. The analysis scripts in bin/python/ read the
// published results layout (real_intervals/, sensitivity/, simulations/,
// benchmarked_calls/), so each task first rebuilds the part of that layout it
// needs from its staged inputs. `manifest` holds one "<relative path>\t<file
// name>" line per file, in the same order as `files`, which are staged one per
// numbered directory because many of them share a name (params.json, and the
// same bench prefix under different settings).

// Environment for every post-hoc task: no bytecode written into bin/python (a new
// entry there would change Manta's cache key), and numerical libraries held to the
// task's CPUs. Unset, they start one thread per core of the node, which on a
// one-CPU allocation slowed the fidelity step roughly sixfold.
def thread_limits(cpus) {
    """
    export PYTHONDONTWRITEBYTECODE=1
    export OMP_NUM_THREADS=${cpus} OPENBLAS_NUM_THREADS=${cpus} MKL_NUM_THREADS=${cpus} NUMEXPR_NUM_THREADS=${cpus}
    """
}

// Shell code that links every staged file to its place under results/.
def link_results_tree(manifest, cpus) {
    """
    ${thread_limits(cpus)}
    i=0
    while IFS=\$'\\t' read -r rel name; do
        [ -n "\$rel" ] || continue
        i=\$((i + 1))
        src="f\$i/\$name"
        [ -e "\$src" ] || src="f/\$name"
        mkdir -p "results/\$(dirname "\$rel")"
        ln -s "\$PWD/\$src" "results/\$rel"
    done <<'MANIFEST'
${manifest}
MANIFEST
    """
}

process POSTHOC_METRICS {
    tag "${assembly}"
    label 'process_single'

    input:
    tuple val(manifest), path(files, stageAs: 'f*/*')
    val assembly
    path compare_results

    output:
    path "metrics.tsv"           , emit: metrics
    path "truvari_parameters.tsv", emit: parameters
    path "compare*"              , emit: comparison, optional: true

    script:
    """
    ${link_results_tree(manifest, task.cpus)}
    python3 ${projectDir}/bin/python/collect_benchmark_metrics.py --results results --assembly ${assembly} \\
        --include-simulations --output metrics.tsv --params-output truvari_parameters.tsv
    if [ -n "${compare_results}" ]; then
        python3 ${projectDir}/bin/python/collect_benchmark_metrics.py --results ${compare_results} \\
            --assembly ${assembly} --include-simulations --output compare_baseline_metrics.tsv
        python3 ${projectDir}/bin/python/compare_runs.py --old compare_baseline_metrics.tsv \\
            --new metrics.tsv --prefix compare
    fi
    """

    stub:
    """
    touch metrics.tsv truvari_parameters.tsv
    """
}

process POSTHOC_SVTYPE_ACCOUNTING {
    tag "${assembly}"
    label 'process_single'

    input:
    tuple val(manifest), path(files, stageAs: 'f*/*')
    path hci_bed
    tuple path(truth_vcf), path(truth_tbi)
    val assembly

    output:
    path "svtype_accounting.tsv", emit: table
    path "svtype_accounting.log", emit: log

    script:
    """
    ${link_results_tree(manifest, task.cpus)}
    python3 ${projectDir}/bin/python/svtype_accounting.py --results results --assembly ${assembly} \\
        --hci-bed ${hci_bed} --truth-vcf ${truth_vcf} --output svtype_accounting.tsv \\
        | tee svtype_accounting.log
    """

    stub:
    """
    touch svtype_accounting.tsv svtype_accounting.log
    """
}

process POSTHOC_DECOMPOSITION {
    tag "${pipeline}"
    label 'process_single'

    input:
    tuple val(pipeline), val(manifest), path(files, stageAs: 'f*/*')
    tuple path(target_bed), path(gene_panel_bed)
    val assembly

    output:
    path "*.decomposition.tsv"        , emit: decomposition
    path "*.strata.tsv"               , emit: strata
    path "*.candidate_transitions.tsv", emit: candidate_transitions
    path "*.uncertainty.tsv"          , emit: uncertainty

    script:
    def name = pipeline.replace(':', '_')
    """
    ${link_results_tree(manifest, task.cpus)}
    python3 ${projectDir}/bin/python/metric_decomposition.py --results results --assembly ${assembly} \\
        --pipeline "${pipeline}" --hci-target ${params.transition_hci_target} \\
        --target "${params.transition_target}=${target_bed}" --target "gene_panel=${gene_panel_bed}" \\
        --simulations --prefix ${name}
    python3 ${projectDir}/bin/python/bootstrap_metrics.py --results results --assembly ${assembly} \\
        --pipeline "${pipeline}" --target ${params.transition_target} --target-bed ${target_bed} \\
        --output ${name}.uncertainty.tsv
    """

    stub:
    def name = pipeline.replace(':', '_')
    """
    touch ${name}.decomposition.tsv ${name}.strata.tsv ${name}.candidate_transitions.tsv ${name}.uncertainty.tsv
    """
}

process POSTHOC_STRATIFY {
    tag "${pipeline}"
    label 'process_single'

    input:
    tuple val(pipeline), val(manifest), path(files, stageAs: 'f*/*')
    tuple path(target_bed), path(gene_panel_bed)
    val assembly

    output:
    path "*.stratified.tsv", emit: stratified

    script:
    def (technology, caller) = pipeline.tokenize(':')
    def name = pipeline.replace(':', '_')
    """
    ${link_results_tree(manifest, task.cpus)}
    python3 ${projectDir}/bin/python/stratify_post_matching.py \\
        --hci-bench results/real_intervals/${technology}/truvari/${caller}/${params.transition_hci_target} \\
        --target "${params.transition_target}=${target_bed}" --target "gene_panel=${gene_panel_bed}" \\
        --assembly ${assembly} --pipeline "${technology} ${caller}" --output ${name}.stratified.tsv
    """

    stub:
    """
    touch ${pipeline.replace(':', '_')}.stratified.tsv
    """
}

// Post-matching stratification and the decomposition's composition-only terms
// score the same truth records and the same HCI false positives, so they must
// agree exactly. They may differ in one way only: a 30-49 bp candidate that
// matched in HCI stays a TP after stratification, while an independently
// restricted benchmark drops it (neither TP nor FP) once its truth partner is
// outside the target. Stratified TP-comp may therefore exceed the
// decomposition's HCI-TP candidates, and never fall below them.
process POSTHOC_STRATIFICATION_CHECK {
    label 'process_single'

    input:
    path decompositions
    path stratified

    output:
    path "stratification_check.txt", emit: report

    script:
    """
    python3 - <<'EOF' | tee stratification_check.txt
import csv, sys
from pathlib import Path
bad = 0
for strat in sorted(Path('.').glob('*.stratified.tsv')):
    name = strat.name[:-len('.stratified.tsv')]
    dec = Path(f'{name}.decomposition.tsv')
    comp = {r['target']: r for r in csv.DictReader(dec.open(), delimiter='\\t')}
    for r in csv.DictReader(strat.open(), delimiter='\\t'):
        d = comp[r['target']]
        extra_tp = int(r['TP-comp']) - int(d['candidate_hci_tp'])
        ok = (int(r['TP-base']) == int(d['truth_hci_tp'])
              and int(r['truth_denominator']) == int(d['n_truth'])
              and int(r['FP']) == int(d['candidate_hci_fp'])
              and extra_tp >= 0)
        bad += not ok
        print(f"{r['pipeline']} {r['target']}: {'ok' if ok else 'MISMATCH'}"
              f" (match-only candidates kept as TP by stratification: {extra_tp})")
sys.exit(1 if bad else 0)
EOF
    """

    stub:
    """
    touch stratification_check.txt
    """
}

process POSTHOC_SENSITIVITY {
    tag "${assembly}"
    label 'process_single'

    input:
    tuple val(manifest), path(files, stageAs: 'f*/*')
    path primary_transitions
    path target_bed
    val threshold_settings
    val recovery_settings
    val assembly

    output:
    path "${assembly}.*.tsv", emit: tables

    script:
    """
    ${link_results_tree(manifest, task.cpus)}
    python3 ${projectDir}/bin/python/sensitivity_transitions.py --results results --assembly ${assembly} \\
        --target ${params.transition_target} --hci-target ${params.transition_hci_target} \\
        --target-label '${params.transition_target_label}' --target-bed ${target_bed} \\
        --threshold-settings '${threshold_settings}' --recovery-settings '${recovery_settings}' \\
        --primary-transitions ${primary_transitions} --prefix ${assembly}
    """

    stub:
    """
    touch ${assembly}.recovery.tsv
    """
}

process POSTHOC_COMPOSITION {
    label 'process_single'

    input:
    path strata

    output:
    path "composition_standardisation.tsv", emit: table

    script:
    // A single staged file arrives as a Path, which Groovy would iterate by name element.
    def strata_args = [strata].flatten().collect { file -> "--strata ${file}" }.join(' ')
    """
    ${thread_limits(task.cpus)}
    python3 ${projectDir}/bin/python/composition_standardisation.py ${strata_args} \\
        --target ${params.transition_target} --output composition_standardisation.tsv
    """

    stub:
    """
    touch composition_standardisation.tsv
    """
}

process POSTHOC_FIDELITY {
    label 'process_single'

    input:
    path target_bed
    path simulated_beds, stageAs: 'simulated_targets/*'
    tuple path(fasta), path(fai)
    tuple val(labels), path(annotations, stageAs: 'annotations/*')

    output:
    path "fidelity.*.tsv", emit: tables

    script:
    // labels and annotations are parallel lists; the label names the column.
    def annotation_args = [[labels].flatten(), [annotations].flatten()].transpose()
        .collect { label, bed -> "--annotation ${label}=${bed}" }.join(' ')
    """
    ${thread_limits(task.cpus)}
    python3 ${projectDir}/bin/python/simulation_fidelity.py --target-bed ${target_bed} \\
        --simulation-dir simulated_targets --reference ${fasta} ${annotation_args} --prefix fidelity
    """

    stub:
    """
    touch fidelity.summary.tsv
    """
}
