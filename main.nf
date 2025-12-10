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
          
        Simulation Options:
          --simulate_targets     Enable target region simulation (default: false)
          --num_simulations      Number of simulations to run (default: 100)
          --gencode_gtf          GENCODE GTF annotation file for simulation
          
        Analysis Options:
          --gather_statistics    Generate statistics and plots (default: false)
        
        Profiles:
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
    if (params.illumina_wes_bam) technologies << "Illumina WES (Manta)"
    if (params.illumina_wgs_bam) technologies << "Illumina WGS (Manta)"
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
    // SUBWORKFLOW: Benchmarking
    //
    ch_truvari_results = Channel.empty()
    if (params.benchmark_vcf && !params.skip_benchmarking) {
        BENCHMARKING(
            SV_CALLING.out.vcfs,
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
    // SUBWORKFLOW: Simulation and benchmarking (optional)
    //
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
            SV_CALLING.out.vcfs,
            params.num_simulations,
            ch_wes_utr,
            ch_high_confidence
        )
        ch_truvari_results = ch_truvari_results.mix(SIMULATE_AND_BENCHMARK.out.truvari_results)
        
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
