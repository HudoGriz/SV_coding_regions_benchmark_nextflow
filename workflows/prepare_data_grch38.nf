#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
    GRCh38 Data Preparation Workflow
========================================================================================
    Downloads and prepares reference data and input files for GRCh38
========================================================================================
*/

// Import nf-core modules
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_EXOME } from '../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_PANEL } from '../modules/nf-core/bedtools/intersect/main'

// Import local modules
include { DOWNLOAD_ILLUMINA_WGS_BAM_GRCH38 } from '../modules/local/download_bam'
include { DOWNLOAD_PACBIO_BAM_GRCH38 } from '../modules/local/download_bam'
include { DOWNLOAD_ONT_BAM_GRCH38 } from '../modules/local/download_bam'
include { DOWNLOAD_REFERENCE_GENOME_GRCH38 } from '../modules/local/download_reference'
include { DOWNLOAD_GIAB_TRUTH_SET_GRCH38 } from '../modules/local/download_truth_set'
include { DOWNLOAD_GIAB_TRUTH_SET_GRCH37_LIFTOVER } from '../modules/local/download_truth_set'
include { DOWNLOAD_TANDEM_REPEATS_GRCH38 } from '../modules/local/download_annotations'
include { DOWNLOAD_GENCODE_GTF_GRCH38 } from '../modules/local/download_annotations'
include { CREATE_EXOME_UTR_BED } from '../modules/local/create_target_beds'

workflow PREPARE_DATA_GRCH38 {
    
    take:
    project_dir    // Path to project directory
    
    main:
    
    // Define output directories
    data_dir = "${project_dir}/data"
    references_dir = "${data_dir}/references"
    
    //
    // Download BAM files (GRCh38)
    //
    
    // Illumina WGS
    ch_illumina_wgs_urls = Channel.value([
        [id: 'Illumina_WGS_GRCh38', technology: 'Illumina_WGS'],
        'https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam',
        "${data_dir}/Illumina_wgs/bam_GRCh38"
    ])
    DOWNLOAD_ILLUMINA_WGS_BAM_GRCH38(ch_illumina_wgs_urls)
    
    // PacBio
    ch_pacbio_urls = Channel.value([
        [id: 'PacBio_GRCh38', technology: 'PacBio'],
        'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam',
        "${data_dir}/Pacbio/bam_GRCh38"
    ])
    DOWNLOAD_PACBIO_BAM_GRCH38(ch_pacbio_urls)
    
    // ONT
    ch_ont_urls = Channel.value([
        [id: 'ONT_GRCh38', technology: 'ONT'],
        'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/HG002_GRCh38_ONT-UL_UCSC_20200508.phased.bam',
        "${data_dir}/ONT/bam_GRCh38"
    ])
    DOWNLOAD_ONT_BAM_GRCH38(ch_ont_urls)
    
    //
    // Download reference files
    //
    
    // Reference genome (GRCh38 no alt)
    ch_ref_url = Channel.value([
        [id: 'GRCh38_no_alt'],
        'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz',
        references_dir
    ])
    DOWNLOAD_REFERENCE_GENOME_GRCH38(ch_ref_url)
    
    // Index reference genome
    SAMTOOLS_FAIDX(
        DOWNLOAD_REFERENCE_GENOME_GRCH38.out.fasta,
        [[],[]]
    )
    
    // GIAB SV Truth Set (GRCh38)
    ch_truth_grch38_urls = Channel.value([
        [id: 'HG002_T2TQ100_V1.0_stvar_GRCh38'],
        'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.012-20231107',
        references_dir,
        'GRCh38'
    ])
    DOWNLOAD_GIAB_TRUTH_SET_GRCH38(ch_truth_grch38_urls)
    
    // GIAB SV Truth Set (GRCh37 liftover)
    ch_truth_grch37_urls = Channel.value([
        [id: 'HG002_T2TQ100_V1.0_stvar_GRCh37'],
        'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.012-20231107',
        references_dir,
        'GRCh37'
    ])
    DOWNLOAD_GIAB_TRUTH_SET_GRCH37_LIFTOVER(ch_truth_grch37_urls)
    
    // Tandem Repeats (GRCh38)
    ch_trf_url = Channel.value([
        [id: 'GRCh38_trf'],
        'https://raw.githubusercontent.com/PacificBiosciences/pbsv/refs/heads/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed',
        references_dir
    ])
    DOWNLOAD_TANDEM_REPEATS_GRCH38(ch_trf_url)
    
    // GENCODE GTF (latest release for GRCh38)
    ch_gencode_url = Channel.value([
        [id: 'gencode_v49'],
        'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v49.annotation.gtf.gz',
        references_dir
    ])
    DOWNLOAD_GENCODE_GTF_GRCH38(ch_gencode_url)
    
    //
    // Create target BED files
    //
    
    // Create exome + UTR BED file using R script
    CREATE_EXOME_UTR_BED(
        DOWNLOAD_GENCODE_GTF_GRCH38.out.gtf,
        Channel.value(file("${project_dir}/singularity_images/r-env_4-4-1.sif")),
        'GRCh38'
    )
    
    //
    // Intersect target sets with truth set
    //
    
    // Intersect exome+UTR with GRCh38 truth set
    // nf-core BEDTOOLS_INTERSECT takes two inputs:
    // 1. tuple val(meta), path(intervals1), path(intervals2)
    // 2. tuple val(meta2), path(chrom_sizes) - empty for BED files
    ch_exome_intersect = DOWNLOAD_GIAB_TRUTH_SET_GRCH38.out.bed
        .combine(CREATE_EXOME_UTR_BED.out.bed)
        .map { bed1, bed2 -> 
            [[id: 'exome_utr_GRCh38_HG002_T2TQ100'], bed1, bed2]
        }
    
    BEDTOOLS_INTERSECT_EXOME(
        ch_exome_intersect,
        [[id: 'null'], []]  // Empty chrom_sizes channel
    )
    
    // Intersect Paediatric disorders with GRCh38 truth set
    ch_paediatric_bed = Channel.value(file("${references_dir}/Paediatric_disorders_GRCh38.bed"))
    
    ch_panel_intersect = DOWNLOAD_GIAB_TRUTH_SET_GRCH38.out.bed
        .combine(ch_paediatric_bed)
        .map { bed1, bed2 -> 
            [[id: 'paediatric_disorders_GRCh38_HG002_T2TQ100'], bed1, bed2]
        }
    
    BEDTOOLS_INTERSECT_PANEL(
        ch_panel_intersect,
        [[id: 'null'], []]  // Empty chrom_sizes channel
    )
    
    emit:
    fasta = DOWNLOAD_REFERENCE_GENOME_GRCH38.out.fasta
    fai = SAMTOOLS_FAIDX.out.fai
    truth_vcf_grch38 = DOWNLOAD_GIAB_TRUTH_SET_GRCH38.out.vcf
    truth_bed_grch38 = DOWNLOAD_GIAB_TRUTH_SET_GRCH38.out.bed
    truth_vcf_grch37 = DOWNLOAD_GIAB_TRUTH_SET_GRCH37_LIFTOVER.out.vcf
    truth_bed_grch37 = DOWNLOAD_GIAB_TRUTH_SET_GRCH37_LIFTOVER.out.bed
    exome_utr_targets = BEDTOOLS_INTERSECT_EXOME.out.intersect
    gene_panel_targets = BEDTOOLS_INTERSECT_PANEL.out.intersect
    tandem_repeats = DOWNLOAD_TANDEM_REPEATS_GRCH38.out.bed
    illumina_wgs_bam = DOWNLOAD_ILLUMINA_WGS_BAM_GRCH38.out.bam
    pacbio_bam = DOWNLOAD_PACBIO_BAM_GRCH38.out.bam
    ont_bam = DOWNLOAD_ONT_BAM_GRCH38.out.bam
}

workflow {
    // Default parameters
    params.project_dir = projectDir.parent
    
    // Run preparation
    PREPARE_DATA_GRCH38(params.project_dir)
}

workflow.onComplete {
    log.info """\
        GRCh38 Data Preparation completed!
        ====================================
        Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
        
        Downloaded files are in:
        - BAM files: ${params.project_dir}/data/*/bam_GRCh38/
        - References: ${params.project_dir}/data/references/
        
        You can now run the main SV calling pipeline with GRCh38 data.
        """
        .stripIndent()
}
