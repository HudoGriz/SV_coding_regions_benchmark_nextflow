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
        
        ch_illumina_wes_bam = Channel.value([
            [id: 'Illumina_WES', technology: 'Illumina_WES', tool: 'Manta'],
            file(params.illumina_wes_bam, checkIfExists: !is_remote),
            file("${params.illumina_wes_bam}.bai", checkIfExists: !is_remote),
            [],  // target_bed
            []   // target_bed_tbi
        ])
        
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
        
        ch_illumina_wgs_bam = Channel.value([
            [id: 'Illumina_WGS', technology: 'Illumina_WGS', tool: 'Manta'],
            file(params.illumina_wgs_bam, checkIfExists: !is_remote),
            file("${params.illumina_wgs_bam}.bai", checkIfExists: !is_remote),
            [],  // target_bed
            []   // target_bed_tbi
        ])
        
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
        ch_pacbio_bam = Channel.value([
            [id: 'PacBio', technology: 'PacBio'],
            file(params.pacbio_bam),
            file("${params.pacbio_bam}.bai")
        ])
        
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
        ch_ont_bam = Channel.value([
            [id: 'ONT', technology: 'ONT'],
            file(params.ont_bam),
            file("${params.ont_bam}.bai")
        ])
        
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
