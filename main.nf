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

// Import helper functions
import WorkflowHelp

// Show help message if requested
if (params.help) {
    log.info WorkflowHelp.helpMessage()
    exit 0
}

// Validate parameters
WorkflowHelp.validateParameters(params, log)

// Print pipeline information
WorkflowHelp.logPipelineInfo(params, log)

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
include { PREPARE_GIAB_RESOURCES } from './workflows/prepare_giab_resources'
include { PREPARE_DATA_COMPLETE_GRCH37 } from './workflows/preparation/prepare_data_complete_grch37'
include { PREPARE_DATA_COMPLETE_GRCH38 } from './workflows/preparation/prepare_data_complete_grch38'

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {
    
    //
    // SUBWORKFLOW: Prepare GIAB resources (optional - minimal)
    //
    if (params.prepare_giab_resources && !params.prepare_complete_data) {
        PREPARE_GIAB_RESOURCES()
        
        log.info """
        =====================================================
        GIAB resources prepared. Outputs available at:
        - Truth VCF: ${params.project_dir}/data/HG002_references/
        - Annotations: ${params.project_dir}/data/references/
        =====================================================
        """.stripIndent()
    }
    
    //
    // SUBWORKFLOW: Complete data preparation (optional - includes everything)
    //
    if (params.prepare_complete_data) {
        if (params.genome == 'hs37d5' || params.genome == 'GRCh37') {
            PREPARE_DATA_COMPLETE_GRCH37(
                params.skip_singularity_download,
                params.skip_bam_download,
                params.skip_reference_download
            )
            
            log.info """
            =====================================================
            Complete GRCh37 data preparation finished!
            
            Downloaded:
            ${params.skip_singularity_download ? '  ✓ Singularity containers (skipped)' : '  ✓ Singularity containers'}
            ${params.skip_bam_download ? '  ✓ BAM files (skipped)' : '  ✓ BAM files (Illumina WES, WGS, PacBio, ONT)'}
            ${params.skip_reference_download ? '  ✓ Reference genome (skipped)' : '  ✓ Reference genome (hs37d5)'}
              ✓ GIAB truth sets (v0.6)
              ✓ Annotations (tandem repeats, GENCODE v19)
              ✓ Target BED files (exome+UTR)
            
            Outputs location: ${params.project_dir}/data/
            =====================================================
            """.stripIndent()
            
        } else if (params.genome == 'GRCh38') {
            PREPARE_DATA_COMPLETE_GRCH38(
                params.skip_bam_download,
                params.skip_reference_download,
                params.download_grch37_liftover
            )
            
            log.info """
            =====================================================
            Complete GRCh38 data preparation finished!
            
            Downloaded:
            ${params.skip_bam_download ? '  ✓ BAM files (skipped)' : '  ✓ BAM files (Illumina WGS, PacBio, ONT)'}
            ${params.skip_reference_download ? '  ✓ Reference genome (skipped)' : '  ✓ Reference genome (GRCh38 no alt)'}
              ✓ GIAB truth sets (T2TQ100-V1.0)
            ${params.download_grch37_liftover ? '  ✓ GIAB GRCh37 liftover truth sets' : ''}
              ✓ Annotations (tandem repeats, GENCODE v49)
              ✓ Target BED files (exome+UTR)
            
            Outputs location: ${params.project_dir}/data/
            =====================================================
            """.stripIndent()
            
        } else {
            error "Unsupported genome: ${params.genome}. Use 'hs37d5', 'GRCh37', or 'GRCh38'"
        }
        
        // Exit after preparation if no input BAMs specified
        if (!params.illumina_wes_bam && !params.illumina_wgs_bam && !params.pacbio_bam && !params.ont_bam) {
            log.info "Data preparation complete. Exiting (no input BAMs specified for SV calling)."
            return
        }
    }
    
    //
    // Validate input parameters for SV calling workflow
    //
    
    // Check if running SV calling mode (not just data preparation)
    def running_sv_calling = !params.prepare_giab_resources && !params.prepare_complete_data
    
    // Exit early if no BAMs provided and not in preparation mode
    if (running_sv_calling && !params.illumina_wes_bam && !params.illumina_wgs_bam && !params.pacbio_bam && !params.ont_bam) {
        log.error """
        =====================================================
        ERROR: No input BAM files specified!
        
        Please provide at least one BAM file:
          --illumina_wes_bam <path>  Illumina WES BAM
          --illumina_wgs_bam <path>  Illumina WGS BAM
          --pacbio_bam <path>        PacBio BAM
          --ont_bam <path>           Oxford Nanopore BAM
        
        Or run data preparation mode:
          --prepare_giab_resources   Prepare minimal GIAB resources
          --prepare_complete_data    Prepare complete dataset
        =====================================================
        """.stripIndent()
        
        System.exit(1)
    }
    
    // Log which technologies are being analyzed
    if (running_sv_calling) {
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
    }
    
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
        // Collect run directory
        ch_run_dir = Channel.value(file(params.outdir))
        
        ANALYSIS_AND_PLOTS(
            ch_truvari_results,
            ch_run_dir
        )
        
        log.info """
        =====================================================
        Statistics and plots generated
        =====================================================
        """.stripIndent()
    }
}

/*
========================================================================================
    WORKFLOW COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info """\
        Pipeline completed at: ${workflow.complete}
        Execution status: ${workflow.success ? 'OK' : 'failed'}
        Execution duration: ${workflow.duration}
        """
        .stripIndent()
}
