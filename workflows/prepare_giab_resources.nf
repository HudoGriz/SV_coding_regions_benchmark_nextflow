/*
========================================================================================
    PREPARE GIAB RESOURCES WORKFLOW
========================================================================================
    Downloads and prepares Genome in a Bottle (GIAB) HG002 resources for validation
----------------------------------------------------------------------------------------
*/

include { DOWNLOAD_GIAB_TRUTH_SET              } from '../modules/local/download_truth_set'
include { DOWNLOAD_GIAB_TRUTH_SET_GRCH38       } from '../modules/local/download_truth_set'
include { DOWNLOAD_GIAB_TRUTH_SET_GRCH37_LIFTOVER } from '../modules/local/download_truth_set'
include { DOWNLOAD_TANDEM_REPEATS              } from '../modules/local/download_annotations'
include { DOWNLOAD_TANDEM_REPEATS_GRCH38       } from '../modules/local/download_annotations'
include { DOWNLOAD_GENCODE_GTF                 } from '../modules/local/download_annotations'
include { DOWNLOAD_GENCODE_GTF_GRCH38          } from '../modules/local/download_annotations'
include { CREATE_EXOME_UTR_BED                 } from '../modules/local/create_target_beds'

workflow PREPARE_GIAB_RESOURCES {
    
    main:
    // Define metadata
    def meta = [id: 'HG002']
    
    // Define output directories
    def giab_output_dir = "${params.project_dir}/data/HG002_references"
    def ref_output_dir = "${params.project_dir}/data/references"
    
    // Determine genome-specific URLs and processes
    if (params.genome == 'hs37d5') {
        // GRCh37/hs37d5 resources
        def giab_base_url = 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37'
        def trf_url = 'https://raw.githubusercontent.com/hall-lab/speedseq/master/annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed'
        def gencode_url = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz'
        
        // Download GIAB truth set
        giab_ch = Channel.of([meta, giab_base_url, giab_output_dir])
        DOWNLOAD_GIAB_TRUTH_SET(giab_ch)
        
        // Download tandem repeats
        trf_ch = Channel.of([meta, trf_url, ref_output_dir])
        DOWNLOAD_TANDEM_REPEATS(trf_ch)
        
        // Download GENCODE GTF
        gencode_ch = Channel.of([meta, gencode_url, ref_output_dir])
        DOWNLOAD_GENCODE_GTF(gencode_ch)
        
        // Create exome + UTR BED
        CREATE_EXOME_UTR_BED(
            DOWNLOAD_GENCODE_GTF.out.gtf,
            params.r_container,
            'GRCh37'
        )
        
        // Emit outputs
        giab_vcf = DOWNLOAD_GIAB_TRUTH_SET.out.vcf
        giab_bed = DOWNLOAD_GIAB_TRUTH_SET.out.bed
        tandem_repeats = DOWNLOAD_TANDEM_REPEATS.out.bed
        gencode_gtf = DOWNLOAD_GENCODE_GTF.out.gtf
        exome_utr_bed = CREATE_EXOME_UTR_BED.out.bed
        
    } else if (params.genome == 'GRCh38') {
        // GRCh38 resources
        def giab_base_url = 'https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38'
        def trf_url = 'https://raw.githubusercontent.com/hall-lab/speedseq/master/annotations/ceph18.b38.lumpy.exclude.2014-01-15.bed'
        def gencode_url = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz'
        
        // Download GIAB truth set
        giab_ch = Channel.of([meta, giab_base_url, giab_output_dir, params.genome])
        DOWNLOAD_GIAB_TRUTH_SET_GRCH38(giab_ch)
        
        // Download tandem repeats
        trf_ch = Channel.of([meta, trf_url, ref_output_dir])
        DOWNLOAD_TANDEM_REPEATS_GRCH38(trf_ch)
        
        // Download GENCODE GTF
        gencode_ch = Channel.of([meta, gencode_url, ref_output_dir])
        DOWNLOAD_GENCODE_GTF_GRCH38(gencode_ch)
        
        // Create exome + UTR BED
        CREATE_EXOME_UTR_BED(
            DOWNLOAD_GENCODE_GTF_GRCH38.out.gtf,
            params.r_container,
            'GRCh38'
        )
        
        // Emit outputs
        giab_vcf = DOWNLOAD_GIAB_TRUTH_SET_GRCH38.out.vcf
        giab_bed = DOWNLOAD_GIAB_TRUTH_SET_GRCH38.out.bed
        tandem_repeats = DOWNLOAD_TANDEM_REPEATS_GRCH38.out.bed
        gencode_gtf = DOWNLOAD_GENCODE_GTF_GRCH38.out.gtf
        exome_utr_bed = CREATE_EXOME_UTR_BED.out.bed
        
    } else {
        error "Unsupported genome: ${params.genome}. Must be 'hs37d5' or 'GRCh38'"
    }
    
    emit:
    giab_vcf
    giab_bed
    tandem_repeats
    gencode_gtf
    exome_utr_bed
}
