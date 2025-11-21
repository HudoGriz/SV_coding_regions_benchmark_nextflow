/*
========================================================================================
    COMPLETE DATA PREPARATION WORKFLOW - GRCh37
========================================================================================
    Downloads and prepares all necessary files for GRCh37/hs37d5:
    - Singularity containers
    - BAM files (Illumina WES, Illumina WGS, PacBio, ONT)
    - Reference genome
    - GIAB truth sets
    - Annotations (tandem repeats, GENCODE)
    - Target BED files
----------------------------------------------------------------------------------------
*/

include { DOWNLOAD_SINGULARITY_IMAGES                    } from '../../modules/local/download_singularity'
include { DOWNLOAD_BAM as DOWNLOAD_ILLUMINA_WES          } from '../../modules/local/download_bam'
include { DOWNLOAD_BAM as DOWNLOAD_ILLUMINA_WGS          } from '../../modules/local/download_bam'
include { DOWNLOAD_BAM as DOWNLOAD_PACBIO                } from '../../modules/local/download_bam'
include { DOWNLOAD_BAM as DOWNLOAD_ONT                   } from '../../modules/local/download_bam'
include { DOWNLOAD_REFERENCE_GENOME as DOWNLOAD_REFERENCE } from '../../modules/local/download_reference'
include { GUNZIP                                } from '../../modules/nf-core/gunzip/main'
include { SAMTOOLS_FAIDX                        } from '../../modules/nf-core/samtools/faidx/main'
include { DOWNLOAD_GIAB_TRUTH_SET               } from '../../modules/local/download_truth_set'
include { DOWNLOAD_TANDEM_REPEATS               } from '../../modules/local/download_annotations'
include { DOWNLOAD_GENCODE_GTF                  } from '../../modules/local/download_annotations'
include { CREATE_EXOME_UTR_BED                  } from '../../modules/local/create_target_beds'
include { BEDTOOLS_INTERSECT as INTERSECT_EXOME } from '../../modules/nf-core/bedtools/intersect/main'
include { BEDTOOLS_INTERSECT as INTERSECT_PANEL } from '../../modules/nf-core/bedtools/intersect/main'

workflow PREPARE_DATA_COMPLETE_GRCH37 {
    
    take:
    skip_singularity  // Boolean: skip singularity download if containers already exist
    skip_bams         // Boolean: skip BAM downloads if files already exist
    skip_reference    // Boolean: skip reference download if file already exists
    
    main:
    
    def meta = [id: 'HG002', genome: 'hs37d5']
    def data_dir = "${params.project_dir}/data"
    def references_dir = "${data_dir}/references"
    
    // =====================================================================
    // STEP 1: Download Singularity Containers
    // =====================================================================
    
    if (!skip_singularity) {
        def container_configs = [
            [name: 'manta_latest', uri: 'docker://dceoy/manta:latest'],
            [name: 'samtools_latest', uri: 'docker://quay.io/biocontainers/samtools:1.19--h50ea8bc_0'],
            [name: 'cutesv_latest', uri: 'docker://quay.io/biocontainers/cutesv:2.1.1--pyhdfd78af_0'],
            [name: 'pbsv_latest', uri: 'docker://quay.io/pacbio/pbsv:2.10.0_build1'],
            [name: 'sniffles_latest', uri: 'docker://quay.io/biocontainers/sniffles:2.4--pyhdfd78af_0'],
            [name: 'bedtools_latest', uri: 'docker://quay.io/biocontainers/bedtools:2.31.1--h13024bc_3'],
            [name: 'truvari_modded', uri: 'library://blazv/benchmark-sv/truvari_modded:latest'],
            [name: 'r-env_4-4-1', uri: 'library://blazv/benchmark-sv/r-env:4-4-1']
        ]
        
        ch_containers = Channel.fromList(container_configs)
        DOWNLOAD_SINGULARITY_IMAGES(ch_containers)
    }
    
    // =====================================================================
    // STEP 2: Download BAM Files
    // =====================================================================
    
    if (!skip_bams) {
        // Illumina WES
        ch_illumina_wes = Channel.of([
            [
                id: 'HG002_Illumina_WES',
                technology: 'Illumina_WES',
                output_dir: "${data_dir}/Illumina_wes/bam"
            ],
            'https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam'
        ])
        DOWNLOAD_ILLUMINA_WES(ch_illumina_wes)
        
        // Illumina WGS
        ch_illumina_wgs = Channel.of([
            [
                id: 'HG002_Illumina_WGS',
                technology: 'Illumina_WGS',
                output_dir: "${data_dir}/Illumina_wgs/bam"
            ],
            'https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.60x.1.bam'
        ])
        DOWNLOAD_ILLUMINA_WGS(ch_illumina_wgs)
        
        // PacBio
        ch_pacbio = Channel.of([
            [
                id: 'HG002_PacBio',
                technology: 'PacBio',
                output_dir: "${data_dir}/Pacbio/bam"
            ],
            'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_HiFi-Revio_20231031/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh37.bam'
        ])
        DOWNLOAD_PACBIO(ch_pacbio)
        
        // ONT
        ch_ont = Channel.of([
            [
                id: 'HG002_ONT',
                technology: 'ONT',
                output_dir: "${data_dir}/ONT/bam"
            ],
            'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/UCSC_Ultralong_OxfordNanopore_Promethion/HG002_GRCh37_ONT-UL_UCSC_20200508.phased.bam'
        ])
        DOWNLOAD_ONT(ch_ont)
    }
    
    // =====================================================================
    // STEP 3: Download and Prepare Reference Genome
    // =====================================================================
    
    if (!skip_reference) {
        ch_reference = Channel.of([
            [id: 'hs37d5', genome: 'hs37d5'],
            'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz',
            references_dir,
            'human_hs37d5.fasta.gz'
        ])
        DOWNLOAD_REFERENCE(ch_reference)
        
        // Gunzip reference
        // nf-core GUNZIP requires tuple [val(meta), path(archive)]
        ch_gunzip_ref = DOWNLOAD_REFERENCE.out.file
            .map { fasta -> [[id: 'hs37d5'], fasta] }
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
    // STEP 4: Download GIAB Truth Set
    // =====================================================================
    
    def giab_base_url = 'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6'
    def giab_output_dir = references_dir
    
    ch_giab = Channel.of([meta, giab_base_url, giab_output_dir])
    DOWNLOAD_GIAB_TRUTH_SET(ch_giab)
    
    // =====================================================================
    // STEP 5: Download Annotations
    // =====================================================================
    
    // Tandem repeats
    def trf_url = 'https://raw.githubusercontent.com/PacificBiosciences/pbsv/master/annotations/human_hs37d5.trf.bed'
    ch_trf = Channel.of([meta, trf_url, references_dir])
    DOWNLOAD_TANDEM_REPEATS(ch_trf)
    
    // GENCODE GTF
    def gencode_url = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz'
    ch_gencode = Channel.of([meta, gencode_url, references_dir])
    DOWNLOAD_GENCODE_GTF(ch_gencode)
    
    // =====================================================================
    // STEP 6: Create Exome+UTR Target BED
    // =====================================================================
    
    CREATE_EXOME_UTR_BED(
        DOWNLOAD_GENCODE_GTF.out.gtf,
        params.r_container,
        'GRCh37'
    )
    
    // =====================================================================
    // STEP 7: Create Intersected BED Files
    // =====================================================================
    
    // Intersect exome+UTR with GIAB truth set
    ch_intersect_exome = DOWNLOAD_GIAB_TRUTH_SET.out.bed
        .combine(CREATE_EXOME_UTR_BED.out.bed)
        .map { giab_bed, exome_bed ->
            [
                [
                    id: 'exome_utr_HG002_SVs',
                    output_name: 'exome_utr_gtf.HG002_SVs_Tier1.bed'
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
    
    // Intersect Paediatric_disorders with GIAB truth set (if exists)
    if (params.paediatric_disorders_bed) {
        ch_intersect_panel = DOWNLOAD_GIAB_TRUTH_SET.out.bed
            .combine(Channel.value(file(params.paediatric_disorders_bed)))
            .map { giab_bed, panel_bed ->
                [
                    [
                        id: 'paediatric_disorders_HG002_SVs',
                        output_name: 'Paediatric_disorders.HG002_SVs_Tier1.bed'
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
    illumina_wes_bam = skip_bams ? Channel.empty() : DOWNLOAD_ILLUMINA_WES.out.bam
    illumina_wgs_bam = skip_bams ? Channel.empty() : DOWNLOAD_ILLUMINA_WGS.out.bam
    pacbio_bam = skip_bams ? Channel.empty() : DOWNLOAD_PACBIO.out.bam
    ont_bam = skip_bams ? Channel.empty() : DOWNLOAD_ONT.out.bam
    
    // Truth sets
    giab_vcf = DOWNLOAD_GIAB_TRUTH_SET.out.vcf
    giab_bed = DOWNLOAD_GIAB_TRUTH_SET.out.bed
    
    // Annotations
    tandem_repeats = DOWNLOAD_TANDEM_REPEATS.out.bed
    gencode_gtf = DOWNLOAD_GENCODE_GTF.out.gtf
    
    // Target BEDs
    exome_utr_bed = CREATE_EXOME_UTR_BED.out.bed
    exome_utr_intersect = INTERSECT_EXOME.out.intersect
    paediatric_intersect = params.paediatric_disorders_bed ? INTERSECT_PANEL.out.intersect : Channel.empty()
}
