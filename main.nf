#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
    SV Calling and Benchmarking Pipeline
========================================================================================
    Pipeline for calling structural variants across multiple sequencing technologies
    and benchmarking results with Truvari
----------------------------------------------------------------------------------------
*/

/*
========================================================================================
    IMPORT MODULES AND WORKFLOWS
========================================================================================
*/

// Sub-workflows
include { PREPARE_REFERENCES } from './workflows/prepare_references'
include { SV_CALLING } from './workflows/sv_calling'
include { BENCHMARKING } from './workflows/benchmarking'
include { SIMULATE_AND_BENCHMARK } from './workflows/simulate_and_benchmark'
include { ANALYSIS_AND_PLOTS } from './workflows/analysis_and_plots'
include { TARGET_TRANSITION_EVIDENCE } from './workflows/target_transition_evidence'
include { SENSITIVITY_BENCHMARKS } from './workflows/sensitivity_benchmarks'
include { EXCLUDE_HOMREF_CALLS } from './modules/local/exclude_homref_calls'
include { POSTHOC_ANALYSES } from './workflows/posthoc_analyses'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    
    // Show help message if requested
    if (params.help) {
        def helpMessage = """
        =====================================================
        SV CALLING AND BENCHMARKING PIPELINE
        =====================================================
        
        Usage:
          nextflow run main.nf -profile <profile> [options]
        
        Required Arguments:
          --fasta                Reference genome FASTA file
        
        Input BAM Files (at least one required):
          --illumina_wes_bam     Illumina WES BAM file
          --illumina_wgs_bam     Illumina WGS BAM file
          --pacbio_bam           PacBio BAM file
          --ont_bam              Oxford Nanopore BAM file
        
        Benchmarking (optional):
          --benchmark_vcf        Truth VCF for benchmarking
          --skip_benchmarking    Skip Truvari benchmarking (default: false)
          --high_confidence_targets  BED file with high confidence regions
          --gene_panel_targets   BED file with gene panel regions
          --wes_utr_targets      BED file with WES UTR regions
        
        Optional Arguments:
          --outdir               Output directory (default: results)
          --run_name             Run name (default: benchmarking_run)
          --tandem_repeats       Tandem repeats BED file (for Sniffles)
          --skip_delly           Skip Delly on Illumina WGS and WES (default: false)
          --delly_exclude        Delly exclude template (telomeres, centromeres)
          --exclude_homref_calls Drop calls genotyped 0/0 before benchmarking (default: true)
          
        Simulation Options:
          --simulate_targets     Enable target region simulation (default: false)
          --num_simulations      Number of simulations to run (default: 100)
          
        Analysis Options:
          --gather_statistics    Generate statistics and plots (default: false)
          --generate_transition_evidence  Audit target-boundary transitions and create figures
          --sensitivity_benchmarks        Re-score the real targets under alternative settings
          --posthoc_analyses     Decomposition, post-matching stratification, bootstrap,
                                 sensitivity audits and simulation fidelity (default: false)

        Profiles:
          study                  Every analysis of the study: 500 simulated sets, statistics,
                                 transition evidence, sensitivity benchmarks, post-hoc analyses
          test_nfcore            Run with nf-core test data
          test                   Run with minimal test data
          docker                 Use Docker containers
          singularity            Use Singularity containers
        =====================================================
        """.stripIndent()
        
        log.info helpMessage
        return
    }
    
    //
    // Validate input parameters for SV calling workflow
    //
    
    // Exit early if no BAMs provided
    if (!params.illumina_wes_bam && !params.illumina_wgs_bam && !params.pacbio_bam && !params.ont_bam) {
        log.error """
        =====================================================
        ERROR: No input BAM files specified!
        
        Please provide at least one BAM file:
          --illumina_wes_bam <path>  Illumina WES BAM
          --illumina_wgs_bam <path>  Illumina WGS BAM
          --pacbio_bam <path>        PacBio BAM
          --ont_bam <path>           Oxford Nanopore BAM
        =====================================================
        """.stripIndent()
        
        error("No input BAM files specified")
    }
    
    // Log which technologies are being analyzed
    def technologies = []
    if (params.illumina_wes_bam) {
        technologies << (params.skip_delly || !params.wes_sequencing_targets ? "Illumina WES (Manta)" : "Illumina WES (Manta, Delly)")
    }
    if (params.illumina_wgs_bam) {
        technologies << (params.skip_delly ? "Illumina WGS (Manta only - Delly skipped)" : "Illumina WGS (Manta, Delly)")
    }
    if (params.pacbio_bam) {
        if (params.skip_pbsv) {
            technologies << "PacBio (CuteSV only - PBSV skipped)"
        } else {
            technologies << "PacBio (CuteSV, PBSV)"
        }
    }
    if (params.ont_bam) technologies << "ONT (CuteSV, Sniffles)"
    
    log.info """
    =====================================================
    SV CALLING ANALYSIS
    
    Technologies to analyze:
    ${technologies.collect { "  ✓ ${it}" }.join('\n')}
    
    ${!params.illumina_wes_bam && !params.illumina_wgs_bam ? '  ✗ Illumina (no BAM provided - skipping)' : ''}
    ${!params.pacbio_bam ? '  ✗ PacBio (no BAM provided - skipping)' : ''}
    ${!params.ont_bam ? '  ✗ ONT (no BAM provided - skipping)' : ''}
    =====================================================
    """.stripIndent()
    
    //
    // SUBWORKFLOW: Prepare references
    //
    PREPARE_REFERENCES()
    
    ch_fasta = PREPARE_REFERENCES.out.fasta
    ch_fasta_fai = PREPARE_REFERENCES.out.fasta_fai
    ch_benchmark_vcf = PREPARE_REFERENCES.out.benchmark_vcf
    ch_benchmark_vcf_tbi = PREPARE_REFERENCES.out.benchmark_vcf_tbi
    ch_targets = PREPARE_REFERENCES.out.targets
    ch_tandem_repeats = PREPARE_REFERENCES.out.tandem_repeats
    
    //
    // SUBWORKFLOW: SV Calling
    //
    SV_CALLING(
        ch_fasta,
        ch_fasta_fai,
        ch_tandem_repeats
    )

    //
    // Calls genotyped homozygous reference (0/0) are the caller stating that the
    // sample does not carry the variant, and the truth sets count only records
    // that carry an ALT allele, so these calls are dropped before any benchmark.
    // A caller with none keeps its original VCF, so its benchmarks are unchanged
    // and a resumed run reuses them. Only the benchmarking subworkflows read these
    // calls, and all of them need a truth set, so without one the step is skipped.
    //
    ch_calls = SV_CALLING.out.vcfs
    ch_homref_counts = Channel.empty()
    if (params.exclude_homref_calls && params.benchmark_vcf) {
        EXCLUDE_HOMREF_CALLS(SV_CALLING.out.vcfs)
        ch_homref_counts = EXCLUDE_HOMREF_CALLS.out.counts
        ch_calls = SV_CALLING.out.vcfs
            .join(EXCLUDE_HOMREF_CALLS.out.vcf)
            .map { meta, vcf, tbi, filtered_vcf, filtered_tbi, removed ->
                removed.toInteger() > 0 ? [meta, filtered_vcf, filtered_tbi] : [meta, vcf, tbi]
            }
    }
    
    //
    // SUBWORKFLOW: Benchmarking
    //
    ch_truvari_results = Channel.empty()
    if (params.benchmark_vcf && !params.skip_benchmarking) {
        BENCHMARKING(
            ch_calls,
            ch_benchmark_vcf,
            ch_benchmark_vcf_tbi,
            ch_targets,
            ch_fasta,
            ch_fasta_fai
        )
        ch_truvari_results = BENCHMARKING.out.summary
    } else {
        log.info "Skipping Truvari benchmarking (benchmark_vcf=${params.benchmark_vcf}, skip_benchmarking=${params.skip_benchmarking})"
    }
    
    //
    // SUBWORKFLOW: Sensitivity benchmarks on the real targets (optional)
    //
    ch_sensitivity_bench = Channel.empty()
    if (params.sensitivity_benchmarks && params.benchmark_vcf && !params.skip_benchmarking) {
        SENSITIVITY_BENCHMARKS(
            ch_calls,
            ch_targets,
            ch_benchmark_vcf,
            ch_benchmark_vcf_tbi,
            ch_fasta,
            ch_fasta_fai
        )
        ch_sensitivity_bench = SENSITIVITY_BENCHMARKS.out.bench_files
    }

    //
    // SUBWORKFLOW: Simulation and benchmarking (optional)
    //
    ch_simulation_transition_evidence = Channel.empty()
    ch_simulated_evidence_beds = Channel.empty()
    if (params.simulate_targets && params.benchmark_vcf) {
        // Validate required parameters - check all at once
        def missing_params = []
        if (!params.wes_utr_targets) missing_params << "--wes_utr_targets"
        if (!params.high_confidence_targets) missing_params << "--high_confidence_targets"
        
        if (missing_params) {
            error """
            =====================================================
            ERROR: Simulation requires the following parameters:
            ${missing_params.collect { "  ${it} <path/to/file.bed>" }.join('\n')}
            
            Example:
            --wes_utr_targets data/references/exome_utr_gtf.HG002_SVs_Tier1.bed
            --high_confidence_targets data/references/HG002_SVs_Tier1_v0.6.bed
            =====================================================
            """.stripIndent()
        }
        
        // Create channels with file existence validation
        ch_wes_utr = Channel.fromPath(params.wes_utr_targets, checkIfExists: true)
        ch_high_confidence = Channel.fromPath(params.high_confidence_targets, checkIfExists: true)
        
        SIMULATE_AND_BENCHMARK(
            ch_fasta,
            ch_fasta_fai,
            ch_benchmark_vcf,
            ch_benchmark_vcf_tbi,
            ch_calls,
            params.num_simulations,
            ch_wes_utr,
            ch_high_confidence
        )
        ch_truvari_results = ch_truvari_results.mix(SIMULATE_AND_BENCHMARK.out.truvari_results)
        ch_simulation_transition_evidence = SIMULATE_AND_BENCHMARK.out.transition_evidence_input
        ch_simulated_evidence_beds = SIMULATE_AND_BENCHMARK.out.simulated_beds
        
        log.info """
        =====================================================
        Simulating ${params.num_simulations} target sets
        =====================================================
        """.stripIndent()
    }
    
    //
    // SUBWORKFLOW: Analysis and plots (optional)
    //
    if (params.gather_statistics && (params.benchmark_vcf && !params.skip_benchmarking)) {
        ANALYSIS_AND_PLOTS(
            ch_truvari_results
        )
        
        log.info """
        =====================================================
        Statistics and plots generated
        =====================================================
        """.stripIndent()
    }

    if (params.generate_transition_evidence && (params.benchmark_vcf && !params.skip_benchmarking)) {
        TARGET_TRANSITION_EVIDENCE(
            BENCHMARKING.out.transition_evidence_input,
            ch_targets,
            ch_simulation_transition_evidence,
            ch_simulated_evidence_beds
        )

        log.info "Target-transition evidence tables and figures generated"
    }
    
    //
    // SUBWORKFLOW: Post-hoc analyses (optional). They compare the real targets with
    // the simulated interval sets, so they need a truth set and the simulations.
    //
    if (params.posthoc_analyses) {
        if (!params.benchmark_vcf || params.skip_benchmarking || !params.simulate_targets) {
            error "--posthoc_analyses needs --benchmark_vcf, benchmarking enabled and --simulate_targets true"
        }
        if (params.sensitivity_benchmarks && !params.generate_transition_evidence) {
            error "--posthoc_analyses with --sensitivity_benchmarks needs --generate_transition_evidence true"
        }
        POSTHOC_ANALYSES(
            BENCHMARKING.out.bench_files
                .mix(SIMULATE_AND_BENCHMARK.out.bench_files)
                .mix(ch_sensitivity_bench),
            ch_calls,
            SIMULATE_AND_BENCHMARK.out.simulated_beds,
            ch_targets,
            ch_benchmark_vcf.combine(ch_benchmark_vcf_tbi),
            ch_fasta.combine(ch_fasta_fai),
            params.generate_transition_evidence ? TARGET_TRANSITION_EVIDENCE.out.evidence : Channel.empty(),
            params.generate_transition_evidence ? TARGET_TRANSITION_EVIDENCE.out.simulation_evidence : Channel.empty(),
            params.gather_statistics ? ANALYSIS_AND_PLOTS.out.tables : Channel.empty(),
            ch_homref_counts,
            params.reference_assembly
        )
    }

    /*
    ========================================================================================
        WORKFLOW COMPLETION HANDLER
    ========================================================================================
    */
    workflow.onComplete = {
        log.info """
        Pipeline completed at: ${workflow.complete}
        Execution status: ${workflow.success ? 'OK' : 'failed'}
        Execution duration: ${workflow.duration}
        """.stripIndent()
    }
}
