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

include { MANTA as MANTA_WES } from './modules/local/manta'
include { MANTA as MANTA_WGS } from './modules/local/manta'
include { CUTESV as CUTESV_PACBIO } from './modules/local/cutesv'
include { CUTESV as CUTESV_ONT } from './modules/local/cutesv'
include { PBSV } from './modules/local/pbsv'
include { SNIFFLES } from './modules/local/sniffles'
include { TRUVARI_BENCH } from './modules/local/truvari'
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
    // Create input channels
    //
    
    // Reference files
    ch_fasta = Channel.value(file(params.fasta))
    ch_fasta_fai = Channel.value(file("${params.fasta}.fai"))
    
    // Target BED files for benchmarking
    ch_benchmark_vcf = Channel.value(file(params.benchmark_vcf))
    ch_benchmark_vcf_tbi = Channel.value(file("${params.benchmark_vcf}.tbi"))
    
    ch_targets = Channel.from([
        ['high_confidence', file(params.high_confidence_targets)],
        ['gene_panel', file(params.gene_panel_targets)],
        ['wes_utr', file(params.wes_utr_targets)]
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
        ch_illumina_wes_bam = Channel.value([
            [id: 'Illumina_WES', technology: 'Illumina_WES', tool: 'Manta'],
            file(params.illumina_wes_bam),
            file("${params.illumina_wes_bam}.bai")
        ])
        
        MANTA_WES(
            ch_illumina_wes_bam,
            ch_fasta,
            ch_fasta_fai
        )
    }
    
    // Illumina WGS
    if (params.illumina_wgs_bam) {
        ch_illumina_wgs_bam = Channel.value([
            [id: 'Illumina_WGS', technology: 'Illumina_WGS', tool: 'Manta'],
            file(params.illumina_wgs_bam),
            file("${params.illumina_wgs_bam}.bai")
        ])
        
        MANTA_WGS(
            ch_illumina_wgs_bam,
            ch_fasta,
            ch_fasta_fai
        )
    }
    
    // PacBio - CuteSV
    if (params.pacbio_bam) {
        ch_pacbio_bam = Channel.value([
            [id: 'PacBio', technology: 'PacBio'],
            file(params.pacbio_bam),
            file("${params.pacbio_bam}.bai")
        ])
        
        CUTESV_PACBIO(
            ch_pacbio_bam.map { meta, bam, bai -> 
                [[id: meta.id, technology: meta.technology, tool: 'CuteSV'], bam, bai]
            },
            ch_fasta,
            ch_fasta_fai
        )
        
        // PacBio - Pbsv
        PBSV(
            ch_pacbio_bam.map { meta, bam, bai -> 
                [[id: meta.id, technology: meta.technology, tool: 'Pbsv'], bam, bai]
            },
            ch_fasta,
            ch_fasta_fai
        )
    }
    
    // ONT - CuteSV
    if (params.ont_bam) {
        ch_ont_bam = Channel.value([
            [id: 'ONT', technology: 'ONT'],
            file(params.ont_bam),
            file("${params.ont_bam}.bai")
        ])
        
        CUTESV_ONT(
            ch_ont_bam.map { meta, bam, bai -> 
                [[id: meta.id, technology: meta.technology, tool: 'CuteSV'], bam, bai]
            },
            ch_fasta,
            ch_fasta_fai
        )
        
        // ONT - Sniffles
        SNIFFLES(
            ch_ont_bam.map { meta, bam, bai -> 
                [[id: meta.id, technology: meta.technology, tool: 'Sniffles'], bam, bai]
            },
            ch_fasta,
            ch_fasta_fai,
            ch_tandem_repeats
        )
    }
    
    //
    // SUBWORKFLOW: Benchmarking
    //
    
    // Collect all VCF outputs
    ch_all_vcfs = Channel.empty()
    
    if (params.illumina_wes_bam) {
        ch_all_vcfs = ch_all_vcfs.mix(MANTA_WES.out.vcf)
    }
    if (params.illumina_wgs_bam) {
        ch_all_vcfs = ch_all_vcfs.mix(MANTA_WGS.out.vcf)
    }
    if (params.pacbio_bam) {
        ch_all_vcfs = ch_all_vcfs.mix(
            CUTESV_PACBIO.out.vcf,
            PBSV.out.vcf
        )
    }
    if (params.ont_bam) {
        ch_all_vcfs = ch_all_vcfs.mix(
            CUTESV_ONT.out.vcf,
            SNIFFLES.out.vcf
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
    
    // Run Truvari benchmarking
    TRUVARI_BENCH(
        ch_benchmark_input,
        ch_benchmark_vcf,
        ch_benchmark_vcf_tbi,
        ch_fasta,
        ch_fasta_fai
    )
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
