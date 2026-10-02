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

// Precision with and without the candidate inversions every real-target benchmark
// scores as false positives (neither truth set contains inversions).
process POSTHOC_INVERSION_PRECISION {
    tag "${assembly}"
    label 'process_single'

    input:
    tuple val(manifest), path(files, stageAs: 'f*/*')
    path metrics
    val assembly

    output:
    path "inversion_precision.tsv", emit: table

    script:
    """
    ${link_results_tree(manifest, task.cpus)}
    python3 ${projectDir}/bin/python/inversion_precision.py --results results --metrics ${metrics} \\
        --output inversion_precision.tsv
    """

    stub:
    """
    touch inversion_precision.tsv
    """
}

// Truth records admitted by containment and by any overlap: counted directly in
// every simulated set with Truvari's conventions (checked against the real
// target's Truvari counts first), and in the real target by plain VCF-BED
// intersection. The target is clipped to HCI and merged exactly as the
// benchmarks use it.
process POSTHOC_MEMBERSHIP {
    tag "${assembly}"
    label 'process_single'

    input:
    tuple path(tp_base), path(fn)
    tuple path(target_bed), path(hci_bed)
    path simulated_beds, stageAs: 'simulated_targets/*'
    path metrics
    val assembly

    output:
    path "membership.*.tsv", emit: tables

    script:
    """
    ${thread_limits(task.cpus)}
    python3 ${projectDir}/bin/python/pad_target_bed.py --target ${target_bed} --allowed ${hci_bed} \\
        --padding 0 --output target.pad0.bed
    python3 ${projectDir}/bin/python/truth_membership.py --tp-base ${tp_base} --fn ${fn} \\
        --target-bed target.pad0.bed --simulation-dir simulated_targets --metrics ${metrics} \\
        --target ${params.transition_target} --prefix membership
    """

    stub:
    """
    touch membership.summary.tsv membership.simulated.tsv membership.boundary_records.tsv
    """
}

// Second comparator. Truth and candidate records are restricted independently to
// HCI and to the target with Truvari's own membership and filters, so only the
// comparator differs from the Truvari analysis; IDs are rewritten to stable,
// record-derived ones because SVanalyzer identifies variants by ID.
process POSTHOC_SVANALYZER_PREFILTER {
    tag "${name}"
    label 'process_single'

    input:
    tuple val(name), path(vcf), path(tbi)
    tuple path(hci_bed), path(target_bed)

    output:
    tuple val(name), path("${name}.{hci,target}.vcf.gz{,.tbi}"), emit: vcfs

    script:
    def id_prefix = name == 'truth' ? 'truth_' : 'cand_'
    """
    ${thread_limits(task.cpus)}
    python3 ${projectDir}/bin/python/prefilter_vcf.py --vcf ${vcf} --bed ${hci_bed} --overlap 1 \\
        --stable-ids --id-prefix ${id_prefix} --output ${name}.hci.vcf.gz
    python3 ${projectDir}/bin/python/prefilter_vcf.py --vcf ${vcf} --bed ${target_bed} --overlap 1 \\
        --stable-ids --id-prefix ${id_prefix} --output ${name}.target.vcf.gz
    """

    stub:
    """
    touch ${name}.hci.vcf.gz ${name}.hci.vcf.gz.tbi ${name}.target.vcf.gz ${name}.target.vcf.gz.tbi
    """
}

// `svanalyzer benchmark` on HCI and on the target. maxdist mirrors Truvari's
// refdist 500, normsizediff 0.3 its pctsize 0.7 and normdist 1.0 its pctseq 0;
// normshift keeps SVanalyzer's default. SVanalyzer rebuilds <fasta>.fai when it
// is older than the FASTA, so the index is copied, never linked: the copy is newer,
// and a rebuild could never write through a link into the reference directory.
process POSTHOC_SVANALYZER {
    tag "${name}"
    label 'process_single'

    input:
    tuple val(name), path(test_vcfs), path(truth_vcfs)   // <name>.{hci,target}.vcf.gz and truth.*, with indexes
    tuple path(fasta), path(fai)

    output:
    tuple val(name), path("${name}.{hci,target}.{distances,report}"), emit: runs
    path "${name}.{hci,target}.log", emit: logs

    script:
    """
    ln -s ${fasta} reference.fasta
    cp -L ${fai} reference.fasta.fai
    for region in hci target; do
        mkdir -p run_\$region
        ( cd run_\$region && svanalyzer benchmark --ref ../reference.fasta \\
            --test ../${name}.\$region.vcf.gz --truth ../truth.\$region.vcf.gz \\
            --maxdist 500 --normshift 1.0 --normsizediff 0.3 --normdist 1.0 --prefix ../${name}.\$region )
    done
    """

    stub:
    """
    touch ${name}.hci.distances ${name}.hci.report ${name}.target.distances ${name}.target.report \\
        ${name}.hci.log ${name}.target.log
    """
}

// Every HCI-TP to target-FN truth record traced to its HCI partners, per pipeline,
// and one summary table for the assembly in pipeline order. `files` holds the
// SVanalyzer runs and the target-restricted VCFs, truth.target.vcf.gz included.
process POSTHOC_SVANALYZER_ELIGIBILITY {
    tag "${assembly}"
    label 'process_single'

    input:
    val names
    path files, stageAs: 'runs/*'
    val assembly

    output:
    path "svanalyzer_summary.tsv", emit: summary
    path "*.transitions.tsv"     , emit: transitions

    script:
    // names are TECHNOLOGY:CALLER; files carry TECHNOLOGY_CALLER
    def calls = [names].flatten().sort().collect { pipeline ->
        def name = pipeline.replace(':', '_')
        """python3 ${projectDir}/bin/python/svanalyzer_eligibility.py --assembly ${assembly} --pipeline "${pipeline.replace(':', ' ')}" \\
        --hci-prefix runs/${name}.hci --target-prefix runs/${name}.target \\
        --target-truth runs/truth.target.vcf.gz --target-test runs/${name}.target.vcf.gz \\
        --normshift 1.0 --normsizediff 0.3 --normdist 1.0 --summary ${name}.summary.tsv --records ${name}.transitions.tsv"""
    }.join('\n')
    def summaries = [names].flatten().sort().collect { pipeline -> "${pipeline.replace(':', '_')}.summary.tsv" }
    """
    ${thread_limits(task.cpus)}
    ${calls}
    head -n 1 ${summaries[0]} > svanalyzer_summary.tsv
    for f in ${summaries.join(' ')}; do tail -n +2 \$f >> svanalyzer_summary.tsv; done
    """

    stub:
    """
    touch svanalyzer_summary.tsv stub.transitions.tsv
    """
}

// Every value the manuscript and its supplementary tables report, for one
// assembly, from the outputs of the other post-hoc steps and of the statistics,
// transition-evidence and 0/0-removal steps. The inputs are staged into the
// published results layout that manuscript_values.py reads.
process POSTHOC_MANUSCRIPT_VALUES {
    tag "${assembly}"
    label 'process_single'

    input:
    path metrics               , stageAs: 'results/posthoc/*'
    path posthoc_tables        , stageAs: 'results/posthoc/*'
    path decomposition         , stageAs: 'results/posthoc/decomposition/*'
    path stratified            , stageAs: 'results/posthoc/stratified/*'
    path sensitivity           , stageAs: 'results/posthoc/sensitivity/*'
    path membership            , stageAs: 'results/posthoc/membership/*'
    path svanalyzer            , stageAs: 'results/posthoc/svanalyzer/*'
    path statistics            , stageAs: 'results/statistics/tables/*'
    path transitions           , stageAs: 'results/target_transition_evidence/tables/*'
    path simulation_transitions, stageAs: 'results/target_transition_evidence/simulations/tables/*'
    path homref_counts         , stageAs: 'results/benchmarked_calls/*'
    val assembly

    output:
    path "${assembly}.manuscript_numbers.md", emit: numbers
    path "supplementary_tables"            , emit: tables

    script:
    """
    ${thread_limits(task.cpus)}
    python3 ${projectDir}/bin/python/manuscript_values.py --results results --assembly ${assembly} --outdir .
    """

    stub:
    """
    mkdir supplementary_tables
    touch ${assembly}.manuscript_numbers.md
    """
}
