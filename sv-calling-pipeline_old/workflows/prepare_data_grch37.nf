#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
    GRCh37 Data Preparation Workflow
========================================================================================
    Downloads and prepares reference data and input files for GRCh37/hs37d5
========================================================================================
*/

// Import nf-core modules
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_VCF } from '../modules/nf-core/samtools/index/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_EXOME } from '../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_INTERSECT as BEDTOOLS_INTERSECT_PANEL } from '../modules/nf-core/bedtools/intersect/main'

// Import local modules for download tasks
include { DOWNLOAD_SINGULARITY_IMAGES } from '../modules/local/download_singularity'
include { DOWNLOAD_ILLUMINA_WES_BAM } from '../modules/local/download_bam'
include { DOWNLOAD_ILLUMINA_WGS_BAM } from '../modules/local/download_bam'
include { DOWNLOAD_PACBIO_BAM } from '../modules/local/download_bam'
include { DOWNLOAD_ONT_BAM } from '../modules/local/download_bam'
include { DOWNLOAD_REFERENCE_GENOME } from '../modules/local/download_reference'
include { DOWNLOAD_GIAB_TRUTH_SET } from '../modules/local/download_truth_set'
include { DOWNLOAD_TANDEM_REPEATS } from '../modules/local/download_annotations'
include { DOWNLOAD_GENCODE_GTF } from '../modules/local/download_annotations'
include { CREATE_EXOME_UTR_BED } from '../modules/local/create_target_beds'

workflow PREPARE_DATA_GRCH37 {
    
    take:
    project_dir    // Path to project directory
    
    main:
    
    // Define output directories
    data_dir = "${project_dir}/data"
    singularity_dir = "${project_dir}/singularity_images"
    references_dir = "${data_dir}/references"
    
    //
    // Download Singularity images
    //
    DOWNLOAD_SINGULARITY_IMAGES(singularity_dir)
    
    //
    // Download BAM files
    //
    
    // Illumina WES
    ch_illumina_wes_urls = Channel.value([
        [id: 'Illumina_WES', technology: 'Illumina_WES'],
        'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam',
        "${data_dir}/Illumina_wes/bam"
    ])
    DOWNLOAD_ILLUMINA_WES_BAM(ch_illumina_wes_urls)
    
    // Illumina WGS
    ch_illumina_wgs_urls = Channel.value([
        [id: 'Illumina_WGS', technology: 'Illumina_WGS'],
        'https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.60x.1.bam',
        "${data_dir}/Illumina_wgs/bam"
    ])
    DOWNLOAD_ILLUMINA_WGS_BAM(ch_illumina_wgs_urls)
    
    // PacBio
    ch_pacbio_urls = Channel.value([
        [id: 'PacBio', technology: 'PacBio'],
        'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh37.bam',
        "${data_dir}/Pacbio/bam"
    ])
    DOWNLOAD_PACBIO_BAM(ch_pacbio_urls)
    
    // ONT
    ch_ont_urls = Channel.value([
        [id: 'ONT', technology: 'ONT'],
        'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/HG002_GRCh37_ONT-UL_UCSC_20200508.phased.bam',
        "${data_dir}/ONT/bam"
    ])
    DOWNLOAD_ONT_BAM(ch_ont_urls)
    
    //
    // Download reference files
    //
    
    // Reference genome (hs37d5)
    ch_ref_url = Channel.value([
        [id: 'hs37d5'],
        'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz',
        references_dir
    ])
    DOWNLOAD_REFERENCE_GENOME(ch_ref_url)
    
    // Index reference genome
    SAMTOOLS_FAIDX(
        DOWNLOAD_REFERENCE_GENOME.out.fasta,
        [[],[]]  // No FAI or GZI input needed
    )
    
    // GIAB SV Truth Set
    ch_truth_urls = Channel.value([
        [id: 'HG002_SVs_Tier1_v0.6'],
        'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6',
        references_dir
    ])
    DOWNLOAD_GIAB_TRUTH_SET(ch_truth_urls)
    
    // Tandem Repeats
    ch_trf_url = Channel.value([
        [id: 'hs37d5_trf'],
        'https://raw.githubusercontent.com/PacificBiosciences/pbsv/master/annotations/human_hs37d5.trf.bed',
        references_dir
    ])
    DOWNLOAD_TANDEM_REPEATS(ch_trf_url)
    
    // GENCODE GTF
    ch_gencode_url = Channel.value([
        [id: 'gencode_v19'],
        'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz',
        references_dir
    ])
    DOWNLOAD_GENCODE_GTF(ch_gencode_url)
    
    //
    // Create target BED files
    //
    
    // Create exome + UTR BED file using R script
    CREATE_EXOME_UTR_BED(
        DOWNLOAD_GENCODE_GTF.out.gtf,
        DOWNLOAD_SINGULARITY_IMAGES.out.r_env,
        'GRCh37'
    )
    
    //
    // Intersect target sets with truth set
    //
    
    // Intersect exome+UTR with truth set
    BEDTOOLS_INTERSECT_EXOME(
        Channel.value([
            [id: 'exome_utr_HG002_SVs_Tier1'],
            DOWNLOAD_GIAB_TRUTH_SET.out.bed,
            CREATE_EXOME_UTR_BED.out.bed,
            []  // No chromosome sizes needed
        ])
    )
    
    // Intersect Paediatric disorders with truth set
    // Note: Assuming Paediatric_disorders.bed is provided separately
    ch_paediatric_bed = Channel.value(file("${references_dir}/Paediatric_disorders.bed"))
    
    BEDTOOLS_INTERSECT_PANEL(
        Channel.value([
            [id: 'paediatric_disorders_HG002_SVs_Tier1'],
            DOWNLOAD_GIAB_TRUTH_SET.out.bed,
            ch_paediatric_bed,
            []
        ])
    )
    
    emit:
    fasta = DOWNLOAD_REFERENCE_GENOME.out.fasta
    fai = SAMTOOLS_FAIDX.out.fai
    truth_vcf = DOWNLOAD_GIAB_TRUTH_SET.out.vcf
    truth_bed = DOWNLOAD_GIAB_TRUTH_SET.out.bed
    exome_utr_targets = BEDTOOLS_INTERSECT_EXOME.out.bed
    gene_panel_targets = BEDTOOLS_INTERSECT_PANEL.out.bed
    tandem_repeats = DOWNLOAD_TANDEM_REPEATS.out.bed
    illumina_wes_bam = DOWNLOAD_ILLUMINA_WES_BAM.out.bam
    illumina_wgs_bam = DOWNLOAD_ILLUMINA_WGS_BAM.out.bam
    pacbio_bam = DOWNLOAD_PACBIO_BAM.out.bam
    ont_bam = DOWNLOAD_ONT_BAM.out.bam
}

workflow {
    // Default parameters
    params.project_dir = projectDir.parent
    
    // Run preparation
    PREPARE_DATA_GRCH37(params.project_dir)
}

workflow.onComplete {
    log.info """\
        GRCh37 Data Preparation completed!
        ====================================
        Status: ${workflow.success ? 'SUCCESS' : 'FAILED'}
        
        Downloaded files are in:
        - BAM files: ${params.project_dir}/data/*/bam/
        - References: ${params.project_dir}/data/references/
        - Containers: ${params.project_dir}/singularity_images/
        
        You can now run the main SV calling pipeline.
        """
        .stripIndent()
}
