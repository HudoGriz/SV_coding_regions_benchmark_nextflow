/*
========================================================================================
    SV CALLING WORKFLOW
========================================================================================
    Calls structural variants using multiple technologies and tools
----------------------------------------------------------------------------------------
*/

include { MANTA_GERMLINE as MANTA_WES } from '../modules/nf-core/manta/germline/main'
include { MANTA_GERMLINE as MANTA_WGS } from '../modules/nf-core/manta/germline/main'
include { CUTESV as CUTESV_PACBIO } from '../modules/nf-core/cutesv/main'
include { CUTESV as CUTESV_ONT } from '../modules/nf-core/cutesv/main'
include { PBSV_DISCOVER } from '../modules/nf-core/pbsv/discover/main'
include { PBSV_CALL } from '../modules/nf-core/pbsv/call/main'
include { SNIFFLES } from '../modules/nf-core/sniffles/main'
include { TABIX_BGZIPTABIX as BGZIP_TABIX_PBSV } from '../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as BGZIP_TABIX_CUTESV_PACBIO } from '../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as BGZIP_TABIX_CUTESV_ONT } from '../modules/nf-core/tabix/bgziptabix/main'
include { SAMTOOLS_INDEX as INDEX_WES_BAM } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_WGS_BAM } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_PACBIO_BAM } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as INDEX_ONT_BAM } from '../modules/nf-core/samtools/index/main'
include { TABIX_TABIX as TABIX_WES_TARGETS } from '../modules/nf-core/tabix/tabix/main'

workflow SV_CALLING {
    take:
    ch_fasta               // channel: reference FASTA file
    ch_fasta_fai           // channel: reference FAI index
    ch_tandem_repeats      // channel: optional tandem repeats BED (for Sniffles)
    
    main:
    
    // Initialize output channel
    ch_all_vcfs = Channel.empty()
    
    //
    // Illumina WES
    //
    if (params.illumina_wes_bam) {
        def is_remote = params.illumina_wes_bam.startsWith('http://') || 
                       params.illumina_wes_bam.startsWith('https://')
        
        def wes_bam = file(params.illumina_wes_bam, checkIfExists: !is_remote)
        def wes_bai = file("${params.illumina_wes_bam}.bai")
        
        // Check if BAM index exists, create if missing
        if (!wes_bai.exists() && !is_remote) {
            log.info "Creating BAM index for Illumina WES: ${params.illumina_wes_bam}"
            INDEX_WES_BAM(
                Channel.value([[id: 'Illumina_WES'], wes_bam])
            )
            ch_wes_bam_indexed = INDEX_WES_BAM.out.bai.map { meta, bai -> 
                [meta, wes_bam, bai]
            }
        } else {
            ch_wes_bam_indexed = Channel.value([
                [id: 'Illumina_WES'], wes_bam, wes_bai
            ])
        }
        
        // Prepare target BED if provided
        if (params.wes_sequencing_targets) {
            def target_bed = file(params.wes_sequencing_targets, checkIfExists: true)
            def target_tbi = file("${params.wes_sequencing_targets}.tbi")
            
            // Check if tabix index exists, create if missing
            if (!target_tbi.exists()) {
                log.info "Creating tabix index for WES targets: ${params.wes_sequencing_targets}"
                TABIX_WES_TARGETS(
                    Channel.value([[id: 'wes_targets'], target_bed])
                )
                ch_wes_targets = TABIX_WES_TARGETS.out.index.map { meta, index_file -> 
                    [target_bed, index_file]
                }
            } else {
                ch_wes_targets = Channel.value([target_bed, target_tbi])
            }
        } else {
            ch_wes_targets = Channel.value([[], []])
        }
        
        // Combine BAM and targets
        ch_illumina_wes_bam = ch_wes_bam_indexed
            .combine(ch_wes_targets)
            .map { tuple ->
                def meta = tuple[0]
                def bam = tuple[1]
                def bai = tuple[2]
                def target_bed = tuple[3]
                def target_tbi = tuple[4]
                [
                    [id: meta.id, technology: 'Illumina_WES', tool: 'Manta'],
                    bam,
                    bai,
                    target_bed,
                    target_tbi
                ]
            }
        
        MANTA_WES(
            ch_illumina_wes_bam,
            ch_fasta.map { f -> [[id: 'fasta'], f] },
            ch_fasta_fai.map { f -> [[id: 'fai'], f] },
            []  // config file (optional)
        )
        
        ch_all_vcfs = ch_all_vcfs.mix(
            MANTA_WES.out.diploid_sv_vcf.join(MANTA_WES.out.diploid_sv_vcf_tbi)
        )
    }
    
    //
    // Illumina WGS
    //
    if (params.illumina_wgs_bam) {
        def is_remote = params.illumina_wgs_bam.startsWith('http://') || 
                       params.illumina_wgs_bam.startsWith('https://')
        
        def wgs_bam = file(params.illumina_wgs_bam, checkIfExists: !is_remote)
        def wgs_bai = file("${params.illumina_wgs_bam}.bai")
        
        // Check if BAM index exists, create if missing
        if (!wgs_bai.exists() && !is_remote) {
            log.info "Creating BAM index for Illumina WGS: ${params.illumina_wgs_bam}"
            INDEX_WGS_BAM(
                Channel.value([[id: 'Illumina_WGS'], wgs_bam])
            )
            ch_illumina_wgs_bam = INDEX_WGS_BAM.out.bai.map { meta, bai -> 
                [
                    [id: meta.id, technology: 'Illumina_WGS', tool: 'Manta'],
                    wgs_bam,
                    bai,
                    [],  // target_bed
                    []   // target_bed_tbi
                ]
            }
        } else {
            ch_illumina_wgs_bam = Channel.value([
                [id: 'Illumina_WGS', technology: 'Illumina_WGS', tool: 'Manta'],
                wgs_bam,
                wgs_bai,
                [],  // target_bed
                []   // target_bed_tbi
            ])
        }
        
        MANTA_WGS(
            ch_illumina_wgs_bam,
            ch_fasta.map { f -> [[id: 'fasta'], f] },
            ch_fasta_fai.map { f -> [[id: 'fai'], f] },
            []  // config file (optional)
        )
        
        ch_all_vcfs = ch_all_vcfs.mix(
            MANTA_WGS.out.diploid_sv_vcf.join(MANTA_WGS.out.diploid_sv_vcf_tbi)
        )
    }
    
    //
    // PacBio - CuteSV
    //
    if (params.pacbio_bam) {
        def pacbio_bam = file(params.pacbio_bam, checkIfExists: true)
        def pacbio_bai = file("${params.pacbio_bam}.bai")
        
        // Check if BAM index exists, create if missing
        if (!pacbio_bai.exists()) {
            log.info "Creating BAM index for PacBio: ${params.pacbio_bam}"
            INDEX_PACBIO_BAM(
                Channel.value([[id: 'PacBio'], pacbio_bam])
            )
            ch_pacbio_bam = INDEX_PACBIO_BAM.out.bai.map { meta, bai -> 
                [meta + [technology: 'PacBio'], pacbio_bam, bai]
            }
        } else {
            ch_pacbio_bam = Channel.value([
                [id: 'PacBio', technology: 'PacBio'],
                pacbio_bam,
                pacbio_bai
            ])
        }
        
        CUTESV_PACBIO(
            ch_pacbio_bam.map { meta, bam, bai -> 
                [[id: meta.id, technology: meta.technology, tool: 'CuteSV'], bam, bai]
            },
            ch_fasta.map { f -> [[id: 'fasta'], f] }
        )
        
        BGZIP_TABIX_CUTESV_PACBIO(CUTESV_PACBIO.out.vcf)
        
        ch_all_vcfs = ch_all_vcfs.mix(BGZIP_TABIX_CUTESV_PACBIO.out.gz_index)
        
        //
        // PacBio - PBSV (optional, may be skipped for test data)
        //
        if (!params.skip_pbsv) {
            PBSV_DISCOVER(
                ch_pacbio_bam.map { meta, bam, bai -> 
                    [[id: meta.id, technology: meta.technology, tool: 'Pbsv'], bam]
                },
                ch_fasta.map { f -> [[id: 'fasta'], f] }
            )
            
            PBSV_CALL(
                PBSV_DISCOVER.out.svsig,
                ch_fasta.map { f -> [[id: 'fasta'], f] }
            )
            
            BGZIP_TABIX_PBSV(PBSV_CALL.out.vcf)
            
            ch_all_vcfs = ch_all_vcfs.mix(BGZIP_TABIX_PBSV.out.gz_index)
        }
    }
    
    //
    // ONT - CuteSV
    //
    if (params.ont_bam) {
        def ont_bam = file(params.ont_bam, checkIfExists: true)
        def ont_bai = file("${params.ont_bam}.bai")
        
        // Check if BAM index exists, create if missing
        if (!ont_bai.exists()) {
            log.info "Creating BAM index for ONT: ${params.ont_bam}"
            INDEX_ONT_BAM(
                Channel.value([[id: 'ONT'], ont_bam])
            )
            ch_ont_bam = INDEX_ONT_BAM.out.bai.map { meta, bai -> 
                [meta + [technology: 'ONT'], ont_bam, bai]
            }
        } else {
            ch_ont_bam = Channel.value([
                [id: 'ONT', technology: 'ONT'],
                ont_bam,
                ont_bai
            ])
        }
        
        CUTESV_ONT(
            ch_ont_bam.map { meta, bam, bai -> 
                [[id: meta.id, technology: meta.technology, tool: 'CuteSV'], bam, bai]
            },
            ch_fasta.map { f -> [[id: 'fasta'], f] }
        )
        
        BGZIP_TABIX_CUTESV_ONT(CUTESV_ONT.out.vcf)
        
        ch_all_vcfs = ch_all_vcfs.mix(BGZIP_TABIX_CUTESV_ONT.out.gz_index)
        
        //
        // ONT - Sniffles
        //
        SNIFFLES(
            ch_ont_bam.map { meta, bam, bai -> 
                [[id: meta.id, technology: meta.technology, tool: 'Sniffles'], bam, bai]
            },
            ch_fasta.map { f -> [[id: 'fasta'], f] },
            params.tandem_repeats ? 
                ch_tandem_repeats.map { f -> [[id: 'tandem_repeats'], f] } :
                Channel.value([[id: 'null'], []]),
            true,   // vcf_output
            false   // snf_output
        )
        
        ch_all_vcfs = ch_all_vcfs.mix(
            SNIFFLES.out.vcf.join(SNIFFLES.out.tbi)
        )
    }
    
    emit:
    vcfs = ch_all_vcfs  // channel: [meta, vcf, tbi]
}
