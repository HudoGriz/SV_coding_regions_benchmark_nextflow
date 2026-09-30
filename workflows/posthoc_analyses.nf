/*
========================================================================================
    POST-HOC ANALYSES WORKFLOW
========================================================================================
    Analyses of the finished benchmarks, enabled with --posthoc_analyses:

      metrics            every benchmark's metrics and its Truvari parameters in one
                         table; optionally compared row by row with another run
                         (--posthoc_compare_results)
      svtype accounting  caller records by type at each step from the scored VCF to
                         the HCI benchmark
      decomposition      recall and precision change from HCI to each target split
                         into composition and label transitions, with the
                         candidate-side (false-positive origin) audit, per pipeline,
                         plus block-bootstrap intervals and rarefied percentile ranks
      stratification     the HCI benchmark stratified by target after matching, and
                         a check that it equals the composition-only values
      sensitivity        transition audits across the threshold grid, and recovery
                         of the primary losses under --extend and padding
      composition        simulated metrics reweighted to the target's SV mix
      fidelity           how the simulated interval sets differ from the real target

    The scripts in bin/python/ read the published results layout, so each task
    rebuilds the part it needs from its staged inputs (modules/local/posthoc_analyses.nf).
----------------------------------------------------------------------------------------
*/

include {
    POSTHOC_METRICS
    POSTHOC_SVTYPE_ACCOUNTING
    POSTHOC_DECOMPOSITION
    POSTHOC_STRATIFY
    POSTHOC_STRATIFICATION_CHECK
    POSTHOC_SENSITIVITY
    POSTHOC_COMPOSITION
    POSTHOC_FIDELITY
} from '../modules/local/posthoc_analyses'

// Where TRUVARI_BENCH publishes a benchmark, relative to --outdir. Mirrors the
// publishDir closure of TRUVARI_BENCH in conf/modules.config.
def bench_dir(meta) {
    def mode = meta.bnd_mode ?: params.primary_bnd_mode
    def root = mode == params.primary_bnd_mode ? '' : "bnd_sensitivity/${mode}/"
    if (meta.target_set == 'simulated') {
        return "${root}simulations/benchmarks/${meta.technology}/${meta.tool}"
    }
    if (meta.target_set == 'sensitivity') {
        return "${root}sensitivity/${meta.setting}/${meta.technology}/${meta.tool}/${meta.target}"
    }
    return "${root}real_intervals/${meta.technology}/truvari/${meta.tool}/${meta.target}"
}

// [manifest, files] for one task from [relative path, file] pairs, sorted by
// path so that the task inputs, and so its cache key, do not depend on the
// order in which upstream tasks finished.
def as_tree(entries) {
    def sorted = entries.collect { entry -> entry }.sort { a, b -> a[0] <=> b[0] }
    [sorted.collect { entry -> "${entry[0]}\t${entry[1].name}" }.join('\n'), sorted.collect { entry -> entry[1] }]
}

// A benchmark in the primary breakend mode, on a real target, and its pipeline key.
def is_primary(meta) {
    (meta.bnd_mode ?: params.primary_bnd_mode) == params.primary_bnd_mode
}

def is_real(meta) {
    meta.target_set != 'simulated' && meta.target_set != 'sensitivity'
}

def pipeline_of(meta) {
    "${meta.technology}:${meta.tool}".toString()
}

// A comma-separated parameter as a list of trimmed, non-empty strings, exactly as
// the sensitivity workflow reads it, so the setting names below match its own.
def split_values(value) {
    value == null ? [] : value.toString().split(',').collect { v -> v.trim() }.findAll { v -> v }
}

workflow POSTHOC_ANALYSES {
    take:
    ch_bench               // channel: [meta, summary, params.json, tp-base, tbi, tp-comp, tbi, fn, tbi, fp, tbi]
    ch_calls               // channel: [meta, vcf, tbi], the VCFs as scored
    ch_simulated_beds      // channel: simulated interval-set BEDs
    ch_targets             // channel: [target_name, bed]
    ch_truth               // channel: [truth vcf, tbi]
    ch_reference           // channel: [fasta, fai]
    ch_target_transitions  // channel: merged target-transition evidence files
    assembly               // value: GRCh37 or GRCh38

    main:
    def hci = params.transition_hci_target
    def target = params.transition_target

    // Every benchmark file with its path under --outdir: [meta, relative path, file]
    ch_bench_entries = ch_bench.flatMap { row ->
        def meta = row[0]
        def dir = bench_dir(meta)
        def summary = row[1]
        def params_json = row[2]
        def entries = [
            [meta, "${dir}/${summary.name}".toString(), summary],
            [meta, "${dir}/${params_json.parent.name}/params.json".toString(), params_json]
        ]
        row[3..-1].each { file -> entries << [meta, "${dir}/${file.name}".toString(), file] }
        entries
    }

    ch_hci_bed = ch_targets.filter { name, _bed -> name == hci }.map { _name, bed -> bed }.first()
    ch_target_beds = ch_targets.filter { name, _bed -> name == target }.map { _name, bed -> bed }
        .combine(ch_targets.filter { name, _bed -> name == 'gene_panel' }.map { _name, bed -> bed })
        .first()
    ch_target_bed = ch_target_beds.map { target_bed, _gene_panel_bed -> target_bed }

    //
    // Metrics and Truvari parameters of every benchmark
    //
    compare_results = params.posthoc_compare_results ? file(params.posthoc_compare_results, checkIfExists: true) : []
    POSTHOC_METRICS(
        ch_bench_entries
            .filter { _meta, rel, _file -> rel.endsWith('.summary.json') || rel.endsWith('/params.json') }
            .map { _meta, rel, file -> [rel, file] }
            .collect(flat: false)
            .map { entries -> as_tree(entries) },
        assembly,
        compare_results
    )

    //
    // SV-type accounting: the scored VCFs and the primary HCI benchmarks
    //
    ch_call_entries = ch_calls.flatMap { meta, vcf, tbi ->
        def rel = "benchmarked_calls/${meta.technology}/${meta.tool}/${meta.technology}-${meta.tool}.benchmarked.vcf.gz"
        [[rel, vcf], ["${rel}.tbi".toString(), tbi]]
    }
    POSTHOC_SVTYPE_ACCOUNTING(
        ch_bench_entries
            .filter { meta, _rel, _file -> is_real(meta) && is_primary(meta) && meta.target == hci }
            .map { _meta, rel, file -> [rel, file] }
            .mix(ch_call_entries)
            .collect(flat: false)
            .map { entries -> as_tree(entries) },
        ch_hci_bed,
        ch_truth,
        assembly
    )

    //
    // Per pipeline: decomposition and bootstrap (real and simulated benchmarks,
    // simulated BEDs), and post-matching stratification (HCI benchmark only).
    // WES is not scored against the simulated sets and is left out, as there.
    //
    ch_simulated_bed_entries = ch_simulated_beds.flatten()
        .map { bed -> ["simulations/target_regions/${bed.name}".toString(), bed] }
        .collect(flat: false)
        .map { entries -> [entries] }

    ch_pipeline_entries = ch_bench_entries
        .filter { meta, _rel, _file -> is_primary(meta) && meta.target_set != 'sensitivity' && meta.technology != 'Illumina_WES' }
        .map { meta, rel, file -> [pipeline_of(meta), meta, [rel, file]] }

    POSTHOC_DECOMPOSITION(
        ch_pipeline_entries
            .map { pipeline, _meta, entry -> [pipeline, entry] }
            .groupTuple()
            .combine(ch_simulated_bed_entries)
            .map { pipeline, entries, bed_entries -> [pipeline] + as_tree(entries + bed_entries) },
        ch_target_beds,
        assembly
    )

    POSTHOC_STRATIFY(
        ch_pipeline_entries
            .filter { _pipeline, meta, _entry -> is_real(meta) && meta.target == hci }
            .map { pipeline, _meta, entry -> [pipeline, entry] }
            .groupTuple()
            .map { pipeline, entries -> [pipeline] + as_tree(entries) },
        ch_target_beds,
        assembly
    )

    POSTHOC_STRATIFICATION_CHECK(
        POSTHOC_DECOMPOSITION.out.decomposition.collect(),
        POSTHOC_STRATIFY.out.stratified.collect()
    )

    //
    // Sensitivity audits: the sensitivity benchmarks, the primary real-target
    // benchmarks they are compared with, and the primary transition table
    //
    if (params.sensitivity_benchmarks) {
        def threshold_settings = split_values(params.sensitivity_refdist).collect { v -> "refdist${v}" } +
            split_values(params.sensitivity_pctsize).collect { v -> "pctsize${v}" } +
            split_values(params.sensitivity_pctseq).collect { v -> "pctseq${v}" } +
            (params.sensitivity_containment ? ['containment'] : [])
        def recovery_settings = split_values(params.sensitivity_extend).collect { v -> "extend${v}" } +
            split_values(params.sensitivity_padding).collect { v -> "pad${v}" }
        POSTHOC_SENSITIVITY(
            ch_bench_entries
                .filter { meta, _rel, _file -> is_primary(meta) && (meta.target_set == 'sensitivity' || is_real(meta)) }
                .map { _meta, rel, file -> [rel, file] }
                .collect(flat: false)
                .map { entries -> as_tree(entries) },
            ch_target_transitions.flatten().filter { file -> file.name.endsWith('.transitions.tsv') }.first(),
            ch_target_bed,
            threshold_settings.join(','),
            recovery_settings.join(','),
            assembly
        )
    }

    //
    // Composition standardisation of the simulated metrics
    //
    POSTHOC_COMPOSITION(POSTHOC_DECOMPOSITION.out.strata.collect())

    //
    // Simulation fidelity: GIAB v3.3 stratifications for segmental duplications
    // and low mappability (overridable), and the tandem-repeat annotation if given
    //
    def giab = 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3'
    def annotation_labels = ['segdups', 'lowmappability']
    def annotation_files = [
        file(params.posthoc_segdups ?: "${giab}/${assembly}@all/SegmentalDuplications/${assembly}_segdups.bed.gz"),
        file(params.posthoc_lowmappability ?: "${giab}/${assembly}@all/Mappability/${assembly}_lowmappabilityall.bed.gz")
    ]
    if (params.tandem_repeats) {
        annotation_labels << 'tandem_repeats'
        annotation_files << file(params.tandem_repeats, checkIfExists: true)
    }
    POSTHOC_FIDELITY(
        ch_target_bed,
        ch_simulated_beds.flatten().collect(),
        ch_reference,
        channel.value([annotation_labels, annotation_files])
    )

    emit:
    metrics     = POSTHOC_METRICS.out.metrics
    check       = POSTHOC_STRATIFICATION_CHECK.out.report
}
