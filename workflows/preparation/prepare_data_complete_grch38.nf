/*
========================================================================================
    COMPLETE DATA PREPARATION WORKFLOW - GRCh38
========================================================================================
    Downloads and prepares all necessary files for GRCh38:
    - BAM files (Illumina WGS, PacBio, ONT)
    - Reference genome
    - GIAB truth sets (both GRCh38 native and GRCh37 liftover)
    - Annotations (tandem repeats, GENCODE)
    - Target BED files
----------------------------------------------------------------------------------------
*/

include { DOWNLOAD_BAM as DOWNLOAD_ILLUMINA_WGS                    } from '../../modules/local/download_bam'
include { DOWNLOAD_BAM as DOWNLOAD_PACBIO                          } from '../../modules/local/download_bam'
include { DOWNLOAD_BAM as DOWNLOAD_ONT                             } from '../../modules/local/download_bam'
include { DOWNLOAD_REFERENCE_GENOME_GRCH38 as DOWNLOAD_REFERENCE   } from '../../modules/local/download_reference'
include { GUNZIP                                } from '../../modules/nf-core/gunzip/main'
include { SAMTOOLS_FAIDX                        } from '../../modules/nf-core/samtools/faidx/main'
include { DOWNLOAD_GIAB_TRUTH_SET_GRCH38        } from '../../modules/local/download_truth_set'
include { DOWNLOAD_GIAB_TRUTH_SET_GRCH37_LIFTOVER } from '../../modules/local/download_truth_set'
include { DOWNLOAD_TANDEM_REPEATS_GRCH38        } from '../../modules/local/download_annotations'
include { DOWNLOAD_GENCODE_GTF_GRCH38           } from '../../modules/local/download_annotations'
include { CREATE_EXOME_UTR_BED                  } from '../../modules/local/create_target_beds'
include { BEDTOOLS_INTERSECT as INTERSECT_EXOME } from '../../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_INTERSECT as INTERSECT_PANEL } from '../../modules/nf-core/bedtools/intersect/main'

workflow PREPARE_DATA_COMPLETE_GRCH38 {
    
    take:
    skip_bams         // Boolean: skip BAM downloads if files already exist
    skip_reference    // Boolean: skip reference download if file already exists
    download_grch37_liftover // Boolean: also download GRCh37 liftover truth set
    
    main:
    
    def meta = [id: 'HG002', genome: 'GRCh38']
    def data_dir = "${params.project_dir}/data"
    def references_dir = "${data_dir}/references"
    
    // =====================================================================
    // STEP 1: Download BAM Files
    // =====================================================================
    
    if (!skip_bams) {
        // Illumina WGS
        ch_illumina_wgs = Channel.of([
            [
                id: 'HG002_Illumina_WGS_GRCh38',
                technology: 'Illumina_WGS',
                output_dir: "${data_dir}/Illumina_wgs/bam_GRCh38"
            ],
            'https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam'
        ])
        DOWNLOAD_ILLUMINA_WGS(ch_illumina_wgs)
        
        // PacBio
        ch_pacbio = Channel.of([
            [
                id: 'HG002_PacBio_GRCh38',
                technology: 'PacBio',
                output_dir: "${data_dir}/Pacbio/bam_GRCh38"
            ],
            'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam'
        ])
        DOWNLOAD_PACBIO(ch_pacbio)
        
        // ONT
        ch_ont = Channel.of([
            [
                id: 'HG002_ONT_GRCh38',
                technology: 'ONT',
                output_dir: "${data_dir}/ONT/bam_GRCh38"
            ],
            'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/HG002_GRCh38_ONT-UL_UCSC_20200508.phased.bam'
        ])
        DOWNLOAD_ONT(ch_ont)
    }
    
    // =====================================================================
    // STEP 2: Download and Prepare Reference Genome
    // =====================================================================
    
    if (!skip_reference) {
        ch_reference = Channel.of([
            [id: 'GRCh38', genome: 'GRCh38'],
            'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz',
            references_dir,
            'human_GRCh38_no_alt_analysis_set.fasta.gz'
        ])
        DOWNLOAD_REFERENCE(ch_reference)
        
        // Gunzip reference
        // nf-core GUNZIP requires tuple [val(meta), path(archive)]
        ch_gunzip_ref = DOWNLOAD_REFERENCE.out.file
            .map { fasta -> [[id: 'GRCh38'], fasta] }
        GUNZIP(ch_gunzip_ref)
        
        // Index reference
        // nf-core SAMTOOLS_FAIDX requires: tuple [meta, fasta], tuple [meta2, fai], val(get_sizes)
        SAMTOOLS_FAIDX(
            GUNZIP.out.gunzip,
            [[id: 'null'], []],  // Empty fai channel (we're generating it)
            false                // Don't generate .sizes file
        )
    }
    
    // =====================================================================
    // STEP 3: Download GIAB Truth Sets
    // =====================================================================
    
    // GRCh38 native truth set
    def giab_grch38_url = 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.012-20231107'
    ch_giab_grch38 = Channel.of([meta, giab_grch38_url, references_dir, 'GRCh38'])
    DOWNLOAD_GIAB_TRUTH_SET_GRCH38(ch_giab_grch38)
    
    // GRCh37 liftover truth set (optional)
    if (download_grch37_liftover) {
        ch_giab_grch37 = Channel.of([meta, giab_grch38_url, references_dir, 'GRCh37'])
        DOWNLOAD_GIAB_TRUTH_SET_GRCH37_LIFTOVER(ch_giab_grch37)
    }
    
    // =====================================================================
    // STEP 4: Download Annotations
    // =====================================================================
    
    // Tandem repeats
    def trf_url = 'https://raw.githubusercontent.com/PacificBiosciences/pbsv/refs/heads/master/annotations/human_GRCh38_no_alt_analysis_set.trf.bed'
    ch_trf = Channel.of([meta, trf_url, references_dir])
    DOWNLOAD_TANDEM_REPEATS_GRCH38(ch_trf)
    
    // GENCODE GTF (v49 - latest release)
    def gencode_url = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/gencode.v49.annotation.gtf.gz'
    ch_gencode = Channel.of([meta, gencode_url, references_dir])
    DOWNLOAD_GENCODE_GTF_GRCH38(ch_gencode)
    
    // =====================================================================
    // STEP 5: Create Exome+UTR Target BED
    // =====================================================================
    
    CREATE_EXOME_UTR_BED(
        DOWNLOAD_GENCODE_GTF_GRCH38.out.gtf,
        params.r_container,
        'GRCh38'
    )
    
    // =====================================================================
    // STEP 6: Create Intersected BED Files
    // =====================================================================
    
    // Intersect exome+UTR with GIAB GRCh38 benchmark regions
    ch_intersect_exome = DOWNLOAD_GIAB_TRUTH_SET_GRCH38.out.bed
        .combine(CREATE_EXOME_UTR_BED.out.bed)
        .map { giab_bed, exome_bed ->
            [
                [
                    id: 'exome_utr_HG002_SVs_GRCh38',
                    output_name: 'exome_utr_gtf.GRCh38_HG002-T2TQ100-V1.0_stvar.bed'
                ],
                giab_bed,
                exome_bed
            ]
        }
    
    // nf-core BEDTOOLS_INTERSECT requires chrom_sizes input (empty for BED files)
    INTERSECT_EXOME(
        ch_intersect_exome,
        [[id: 'null'], []]  // Empty chrom_sizes channel
    )
    
    // Intersect Paediatric_disorders with GIAB GRCh38 benchmark regions (if exists)
    if (params.paediatric_disorders_bed_grch38) {
        ch_intersect_panel = DOWNLOAD_GIAB_TRUTH_SET_GRCH38.out.bed
            .combine(Channel.value(file(params.paediatric_disorders_bed_grch38)))
            .map { giab_bed, panel_bed ->
                [
                    [
                        id: 'paediatric_disorders_HG002_SVs_GRCh38',
                        output_name: 'Paediatric_disorders_GRCh38.GRCh38_HG002-T2TQ100-V1.0_stvar.bed'
                    ],
                    giab_bed,
                    panel_bed
                ]
            }
        
        INTERSECT_PANEL(
            ch_intersect_panel,
            [[id: 'null'], []]  // Empty chrom_sizes channel
        )
    }
    
    emit:
    // Reference files
    // Extract file from tuple [meta, file]
    reference_fasta = skip_reference ? Channel.empty() : GUNZIP.out.gunzip.map { m, file -> file }
    reference_fai = skip_reference ? Channel.empty() : SAMTOOLS_FAIDX.out.fai.map { m, file -> file }
    
    // BAM files
    illumina_wgs_bam = skip_bams ? Channel.empty() : DOWNLOAD_ILLUMINA_WGS.out.bam
    pacbio_bam = skip_bams ? Channel.empty() : DOWNLOAD_PACBIO.out.bam
    ont_bam = skip_bams ? Channel.empty() : DOWNLOAD_ONT.out.bam
    
    // Truth sets
    giab_vcf_grch38 = DOWNLOAD_GIAB_TRUTH_SET_GRCH38.out.vcf
    giab_bed_grch38 = DOWNLOAD_GIAB_TRUTH_SET_GRCH38.out.bed
    giab_vcf_grch37_liftover = download_grch37_liftover ? DOWNLOAD_GIAB_TRUTH_SET_GRCH37_LIFTOVER.out.vcf : Channel.empty()
    giab_bed_grch37_liftover = download_grch37_liftover ? DOWNLOAD_GIAB_TRUTH_SET_GRCH37_LIFTOVER.out.bed : Channel.empty()
    
    // Annotations
    tandem_repeats = DOWNLOAD_TANDEM_REPEATS_GRCH38.out.bed
    gencode_gtf = DOWNLOAD_GENCODE_GTF_GRCH38.out.gtf
    
    // Target BEDs
    exome_utr_bed = CREATE_EXOME_UTR_BED.out.bed
    exome_utr_intersect = INTERSECT_EXOME.out.intersect
    paediatric_intersect = params.paediatric_disorders_bed_grch38 ? INTERSECT_PANEL.out.intersect : Channel.empty()
}
