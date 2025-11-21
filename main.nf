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

// Help message
def helpMessage() {
    log.info"""
    =====================================================
    SV CALLING AND BENCHMARKING PIPELINE
    =====================================================
    
    Usage:
      nextflow run main.nf -profile <profile> [options]
    
    Required Arguments:
      --fasta                Reference genome FASTA file
      --benchmark_vcf        Truth VCF for benchmarking
      --high_confidence_targets  BED file with high confidence regions
      --gene_panel_targets   BED file with gene panel regions
      --wes_utr_targets      BED file with WES UTR regions
    
    Input BAM Files (at least one required):
      --illumina_wes_bam     Illumina WES BAM file
      --illumina_wgs_bam     Illumina WGS BAM file
      --pacbio_bam           PacBio BAM file
      --ont_bam              Oxford Nanopore BAM file
    
    Optional Arguments:
      --outdir               Output directory (default: results)
      --run_name             Run name (default: benchmarking_run)
      --tandem_repeats       Tandem repeats BED file
      
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
    
    Examples:
      # Run with test data
      nextflow run main.nf -profile test_nfcore,docker
      
      # Run with custom data
      nextflow run main.nf -profile docker \\
        --fasta ref.fa \\
        --illumina_wes_bam sample.bam \\
        --benchmark_vcf truth.vcf.gz \\
        --high_confidence_targets regions.bed
    =====================================================
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

// Print pipeline information
log.info """\
    =====================================================
    SV CALLING AND BENCHMARKING PIPELINE
    =====================================================
    run_name            : ${params.run_name}
    reference           : ${params.fasta}
    outdir              : ${params.outdir}
    =====================================================
    """
    .stripIndent()

/*
========================================================================================
    IMPORT MODULES AND WORKFLOWS
========================================================================================
*/

include { SAMTOOLS_FAIDX } from './modules/nf-core/samtools/faidx/main'
include { MANTA_GERMLINE as MANTA_WES } from './modules/nf-core/manta/germline/main'
include { MANTA_GERMLINE as MANTA_WGS } from './modules/nf-core/manta/germline/main'
include { CUTESV as CUTESV_PACBIO } from './modules/nf-core/cutesv/main'
include { CUTESV as CUTESV_ONT } from './modules/nf-core/cutesv/main'
include { PBSV_DISCOVER } from './modules/nf-core/pbsv/discover/main'
include { PBSV_CALL } from './modules/nf-core/pbsv/call/main'
include { SNIFFLES } from './modules/nf-core/sniffles/main'
include { BGZIP_TABIX as BGZIP_TABIX_PBSV } from './modules/local/bgzip_tabix'
include { BGZIP_TABIX as BGZIP_TABIX_CUTESV_PACBIO } from './modules/local/bgzip_tabix'
include { BGZIP_TABIX as BGZIP_TABIX_CUTESV_ONT } from './modules/local/bgzip_tabix'
include { TRUVARI_BENCH } from './modules/local/truvari'
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
    // Create input channels
    //
    
    // Reference files
    // For remote files (URLs), don't check if exists as they may need to be downloaded
    def is_remote = params.fasta.startsWith('http://') || params.fasta.startsWith('https://') || params.fasta.startsWith('ftp://')
    ch_fasta = Channel.value(file(params.fasta, checkIfExists: !is_remote))
    
    // FAI index - create if needed for remote files or if doesn't exist locally
    def fai_path = "${params.fasta}.fai"
    def fai_file = file(fai_path)
    
    if (is_remote || !fai_file.exists()) {
        // Create FAI index using samtools faidx
        SAMTOOLS_FAIDX(
            ch_fasta.map { f -> [[id: 'reference'], f] },
            [[],[]],  // No existing fai provided, will be created
            false     // get_sizes parameter
        )
        ch_fasta_fai = SAMTOOLS_FAIDX.out.fai.map { meta, fai -> fai }
    } else {
        ch_fasta_fai = Channel.value(fai_file)
    }
    
    // Benchmark VCF and index (optional for testing)
    // For remote files, Nextflow will download them automatically
    if (params.benchmark_vcf && !params.skip_benchmarking) {
        def is_vcf_remote = params.benchmark_vcf.startsWith('http://') || params.benchmark_vcf.startsWith('https://') || params.benchmark_vcf.startsWith('ftp://')
        ch_benchmark_vcf = Channel.value(file(params.benchmark_vcf, checkIfExists: !is_vcf_remote))
        ch_benchmark_vcf_tbi = Channel.value(file("${params.benchmark_vcf}.tbi", checkIfExists: !is_vcf_remote))
    } else {
        ch_benchmark_vcf = Channel.empty()
        ch_benchmark_vcf_tbi = Channel.empty()
    }
    
    // Target BED files - check if remote
    def is_targets_remote = params.high_confidence_targets.startsWith('http://') || params.high_confidence_targets.startsWith('https://') || params.high_confidence_targets.startsWith('ftp://')
    ch_targets = Channel.from([
        ['high_confidence', file(params.high_confidence_targets, checkIfExists: !is_targets_remote)],
        ['gene_panel', file(params.gene_panel_targets, checkIfExists: !is_targets_remote)],
        ['wes_utr', file(params.wes_utr_targets, checkIfExists: !is_targets_remote)]
    ])
    
    // Optional: ONT tandem repeats BED
    ch_tandem_repeats = params.tandem_repeats ? 
        Channel.value(file(params.tandem_repeats)) : 
        Channel.empty()
    
    //
    // SUBWORKFLOW: SV Calling
    //
    
    // Illumina WES
    if (params.illumina_wes_bam) {
        def is_wes_remote = params.illumina_wes_bam.startsWith('http://') || params.illumina_wes_bam.startsWith('https://')
        ch_illumina_wes_bam = Channel.value([
            [id: 'Illumina_WES', technology: 'Illumina_WES', tool: 'Manta'],
            file(params.illumina_wes_bam, checkIfExists: !is_wes_remote),
            file("${params.illumina_wes_bam}.bai", checkIfExists: !is_wes_remote),
            [],  // target_bed
            []   // target_bed_tbi
        ])
        
        // nf-core MANTA_GERMLINE requires: tuple [meta, input, index, target_bed, target_bed_tbi], 
        //                                  tuple [meta2, fasta], tuple [meta3, fai], path(config)
        MANTA_WES(
            ch_illumina_wes_bam,
            ch_fasta.map { f -> [[id: 'fasta'], f] },
            ch_fasta_fai.map { f -> [[id: 'fai'], f] },
            []  // config file (optional)
        )
    }
    
    // Illumina WGS
    if (params.illumina_wgs_bam) {
        def is_wgs_remote = params.illumina_wgs_bam.startsWith('http://') || params.illumina_wgs_bam.startsWith('https://')
        ch_illumina_wgs_bam = Channel.value([
            [id: 'Illumina_WGS', technology: 'Illumina_WGS', tool: 'Manta'],
            file(params.illumina_wgs_bam, checkIfExists: !is_wgs_remote),
            file("${params.illumina_wgs_bam}.bai", checkIfExists: !is_wgs_remote),
            [],  // target_bed
            []   // target_bed_tbi
        ])
        
        MANTA_WGS(
            ch_illumina_wgs_bam,
            ch_fasta.map { f -> [[id: 'fasta'], f] },
            ch_fasta_fai.map { f -> [[id: 'fai'], f] },
            []  // config file (optional)
        )
    }
    
    // PacBio - CuteSV
    if (params.pacbio_bam) {
        ch_pacbio_bam = Channel.value([
            [id: 'PacBio', technology: 'PacBio'],
            file(params.pacbio_bam),
            file("${params.pacbio_bam}.bai")
        ])
        
        // nf-core CUTESV requires: tuple [meta, bam, bai], tuple [meta2, fasta]
        CUTESV_PACBIO(
            ch_pacbio_bam.map { meta, bam, bai -> 
                [[id: meta.id, technology: meta.technology, tool: 'CuteSV'], bam, bai]
            },
            ch_fasta.map { f -> [[id: 'fasta'], f] }
        )
        
        // Compress and index CuteSV output
        BGZIP_TABIX_CUTESV_PACBIO(CUTESV_PACBIO.out.vcf)
        
        // PacBio - Pbsv (requires discover + call)
        // Only run if not skipped (test data may not have proper PacBio headers)
        if (!params.skip_pbsv) {
            // Step 1: Discover SV signatures
            PBSV_DISCOVER(
                ch_pacbio_bam.map { meta, bam, bai -> 
                    [[id: meta.id, technology: meta.technology, tool: 'Pbsv'], bam]
                },
                ch_fasta.map { f -> [[id: 'fasta'], f] }
            )
            
            // Step 2: Call SVs from signatures
            PBSV_CALL(
                PBSV_DISCOVER.out.svsig,
                ch_fasta.map { f -> [[id: 'fasta'], f] }
            )
            
            // Compress and index Pbsv output
            BGZIP_TABIX_PBSV(PBSV_CALL.out.vcf)
        }
    }
    
    // ONT - CuteSV
    if (params.ont_bam) {
        ch_ont_bam = Channel.value([
            [id: 'ONT', technology: 'ONT'],
            file(params.ont_bam),
            file("${params.ont_bam}.bai")
        ])
        
        // nf-core CUTESV requires: tuple [meta, bam, bai], tuple [meta2, fasta]
        CUTESV_ONT(
            ch_ont_bam.map { meta, bam, bai -> 
                [[id: meta.id, technology: meta.technology, tool: 'CuteSV'], bam, bai]
            },
            ch_fasta.map { f -> [[id: 'fasta'], f] }
        )
        
        // Compress and index CuteSV output
        BGZIP_TABIX_CUTESV_ONT(CUTESV_ONT.out.vcf)
        
        // ONT - Sniffles
        // nf-core SNIFFLES requires: tuple [meta, bam, bai], tuple [meta2, fasta], 
        //                            tuple [meta3, tandem_file], vcf_output, snf_output
        SNIFFLES(
            ch_ont_bam.map { meta, bam, bai -> 
                [[id: meta.id, technology: meta.technology, tool: 'Sniffles'], bam, bai]
            },
            ch_fasta.map { f -> [[id: 'fasta'], f] },
            params.tandem_repeats ? 
                ch_tandem_repeats.map { f -> [[id: 'tandem_repeats'], f] } :
                Channel.value([[id: 'null'], []]),
            true,   // vcf_output (output VCF)
            false   // snf_output (don't output SNF)
        )
    }
    
    //
    // SUBWORKFLOW: Benchmarking
    //
    
    // Collect all VCF outputs
    ch_all_vcfs = Channel.empty()
    
    if (params.illumina_wes_bam) {
        // MANTA_GERMLINE outputs diploid_sv_vcf and diploid_sv_vcf_tbi
        ch_all_vcfs = ch_all_vcfs.mix(
            MANTA_WES.out.diploid_sv_vcf
                .join(MANTA_WES.out.diploid_sv_vcf_tbi)
        )
    }
    if (params.illumina_wgs_bam) {
        ch_all_vcfs = ch_all_vcfs.mix(
            MANTA_WGS.out.diploid_sv_vcf
                .join(MANTA_WGS.out.diploid_sv_vcf_tbi)
        )
    }
    if (params.pacbio_bam) {
        // BGZIP_TABIX outputs tuple [meta, vcf, tbi]
        ch_all_vcfs = ch_all_vcfs.mix(
            BGZIP_TABIX_CUTESV_PACBIO.out.vcf
        )
        if (!params.skip_pbsv) {
            ch_all_vcfs = ch_all_vcfs.mix(
                BGZIP_TABIX_PBSV.out.vcf
            )
        }
    }
    if (params.ont_bam) {
        // SNIFFLES outputs vcf and tbi separately
        ch_all_vcfs = ch_all_vcfs.mix(
            BGZIP_TABIX_CUTESV_ONT.out.vcf,
            SNIFFLES.out.vcf.join(SNIFFLES.out.tbi)
        )
    }
    
    // Create combinations of VCFs and target sets for benchmarking
    ch_benchmark_input = ch_all_vcfs
        .combine(ch_targets)
        .map { meta, vcf, vcf_tbi, target_name, target_bed ->
            [
                [
                    id: meta.id,
                    technology: meta.technology,
                    tool: meta.tool,
                    target: target_name
                ],
                vcf,
                vcf_tbi,
                target_bed
            ]
        }
    
    // Run Truvari benchmarking (skip if benchmark_vcf is null or skip_benchmarking is true)
    ch_truvari_results = Channel.empty()
    if (params.benchmark_vcf && !params.skip_benchmarking) {
        TRUVARI_BENCH(
            ch_benchmark_input,
            ch_benchmark_vcf,
            ch_benchmark_vcf_tbi,
            ch_fasta,
            ch_fasta_fai
        )
        ch_truvari_results = TRUVARI_BENCH.out.summary
    } else {
        log.info "Skipping Truvari benchmarking (benchmark_vcf=${params.benchmark_vcf}, skip_benchmarking=${params.skip_benchmarking})"
    }
    
    //
    // SUBWORKFLOW: Simulation and benchmarking (optional)
    //
    if (params.simulate_targets && params.gencode_gtf && params.benchmark_vcf) {
        SIMULATE_AND_BENCHMARK(
            ch_fasta,
            ch_fasta_fai,
            file(params.gencode_gtf, checkIfExists: true),
            ch_benchmark_vcf,
            ch_benchmark_vcf_tbi,
            ch_all_vcfs,
            params.num_simulations
        )
        ch_truvari_results = ch_truvari_results.mix(SIMULATE_AND_BENCHMARK.out.truvari_results)
        
        log.info """
        =====================================================
        Simulation completed: ${params.num_simulations} simulated target sets
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
