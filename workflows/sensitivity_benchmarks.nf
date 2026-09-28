/*
========================================================================================
    SENSITIVITY BENCHMARKS WORKFLOW
========================================================================================
    Re-scores the real targets under alternative benchmark settings, each changing
    one thing relative to the primary configuration:

      refdist<N>, pctsize<X>, pctseq<X>  a matching threshold, on every real target,
                                         so each target can be paired with HCI
                                         scored the same way
      containment                        full-containment membership
                                         (--bench-overlaps 0), on every real target
      extend<N>                          candidate-only boundary extension
                                         (Truvari --extend N) on the boundary target;
                                         the truth set and its denominator do not move
      pad<N>                             the boundary target padded by N bp on both
                                         sides, clipped to HCI; truth and candidates
                                         both move

    Only the primary breakend mode is scored, and WES is left out as it is from the
    empirical null. Results publish under <outdir>/sensitivity/<setting>/ and feed
    no headline table. The audits that compare settings are run afterwards from the
    published files (bin/run_revision_posthoc.sh).
----------------------------------------------------------------------------------------
*/

include { TRUVARI_BENCH } from '../modules/nf-core/truvari/bench/main'
include { PAD_TARGET_BED } from '../modules/local/pad_target_bed'

// A comma-separated parameter as a list of trimmed, non-empty strings. A value
// given on the command line may arrive as a number, hence toString().
def split_values(value) {
    value == null ? [] : value.toString().split(',').collect { it.trim() }.findAll { it }
}

workflow SENSITIVITY_BENCHMARKS {
    take:
    ch_vcfs                 // channel: [meta, vcf, tbi]
    ch_targets              // channel: [target_name, bed]
    ch_benchmark_vcf        // channel: truth VCF
    ch_benchmark_vcf_tbi    // channel: truth VCF index
    ch_fasta                // channel: reference FASTA
    ch_fasta_fai            // channel: reference FAI index

    main:
    def primary = [
        refdist: params.truvari_refdist,
        pctsize: params.truvari_pctsize,
        pctovl : params.truvari_pctovl,
        pctseq : params.truvari_pctseq
    ]
    def all_targets = split_values(params.sensitivity_targets)
    def boundary_target = params.transition_target

    def settings = []
    split_values(params.sensitivity_refdist).each { v ->
        settings << [name: "refdist${v}", targets: all_targets, thresholds: primary + [refdist: v]]
    }
    split_values(params.sensitivity_pctsize).each { v ->
        settings << [name: "pctsize${v}", targets: all_targets, thresholds: primary + [pctsize: v]]
    }
    split_values(params.sensitivity_pctseq).each { v ->
        settings << [name: "pctseq${v}", targets: all_targets, thresholds: primary + [pctseq: v]]
    }
    if (params.sensitivity_containment) {
        settings << [name: 'containment', targets: all_targets, thresholds: primary, bench_overlaps: 0]
    }
    split_values(params.sensitivity_extend).each { v ->
        settings << [name: "extend${v}", targets: [boundary_target], thresholds: primary, extend: v]
    }

    // [target_name, setting, bed] for every setting that uses a real target as is
    ch_setting_beds = Channel
        .fromList(settings.collectMany { setting -> setting.targets.collect { target -> [target, setting] } })
        .combine(ch_targets, by: 0)

    // [target_name, setting, bed] for the padded boundary target
    ch_allowed = ch_targets
        .filter { target_name, bed -> target_name == params.transition_hci_target }
        .map { target_name, bed -> bed }

    PAD_TARGET_BED(
        ch_targets
            .filter { target_name, bed -> target_name == boundary_target }
            .combine(Channel.fromList(split_values(params.sensitivity_padding)))
            .combine(ch_allowed)
    )

    ch_padded_beds = PAD_TARGET_BED.out.bed.map { target_name, padding, bed ->
        [target_name, [name: "pad${padding}", thresholds: primary], bed]
    }

    ch_bench_input = ch_vcfs
        .filter { meta, vcf, tbi -> meta.technology != 'Illumina_WES' }
        .combine(ch_setting_beds.mix(ch_padded_beds))
        .combine(ch_benchmark_vcf)
        .combine(ch_benchmark_vcf_tbi)
        .map { meta, vcf, tbi, target_name, setting, bed, truth_vcf, truth_tbi ->
            def t = setting.thresholds
            def truvari_args = "--refdist ${t.refdist} --pctsize ${t.pctsize} --pctovl ${t.pctovl} --pctseq ${t.pctseq}"
            if (setting.extend != null) {
                truvari_args += " --extend ${setting.extend}"
            }
            def bench_meta = meta + [
                target        : target_name,
                target_set    : 'sensitivity',
                setting       : setting.name,
                bnd_mode      : params.primary_bnd_mode,
                bench_overlaps: setting.bench_overlaps != null ? setting.bench_overlaps : 1,
                truvari_args  : truvari_args
            ]
            bench_meta.id = "${meta.technology}_${meta.tool}_${target_name}_${setting.name}"
            [bench_meta, vcf, tbi, truth_vcf, truth_tbi, bed]
        }

    TRUVARI_BENCH(
        ch_bench_input,
        ch_fasta.map { fasta -> [[id: 'reference'], fasta] },
        ch_fasta_fai.map { fai -> [[id: 'reference'], fai] }
    )

    emit:
    summary     = TRUVARI_BENCH.out.summary  // channel: [meta, summary.json]
    padded_beds = PAD_TARGET_BED.out.bed     // channel: [target_name, padding, bed]
}
