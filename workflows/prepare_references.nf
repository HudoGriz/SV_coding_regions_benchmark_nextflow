/*
========================================================================================
    PREPARE REFERENCES WORKFLOW
========================================================================================
    Prepares reference files and creates necessary indices
----------------------------------------------------------------------------------------
*/

include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'

workflow PREPARE_REFERENCES {
    
    main:
    
    // Reference FASTA file
    if (!params.fasta) {
        error "ERROR: Reference FASTA file must be specified with --fasta parameter"
    }
    
    // For remote files (URLs), don't check if exists as they may need to be downloaded
    def is_remote = params.fasta.startsWith('http://') || 
                    params.fasta.startsWith('https://') || 
                    params.fasta.startsWith('ftp://')
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
        def is_vcf_remote = params.benchmark_vcf.startsWith('http://') || 
                           params.benchmark_vcf.startsWith('https://') || 
                           params.benchmark_vcf.startsWith('ftp://')
        ch_benchmark_vcf = Channel.value(file(params.benchmark_vcf, checkIfExists: !is_vcf_remote))
        ch_benchmark_vcf_tbi = Channel.value(file("${params.benchmark_vcf}.tbi", checkIfExists: !is_vcf_remote))
    } else {
        ch_benchmark_vcf = Channel.empty()
        ch_benchmark_vcf_tbi = Channel.empty()
    }
    
    // Target BED files - check if remote
    if (params.high_confidence_targets && params.gene_panel_targets && params.wes_utr_targets) {
        def is_targets_remote = params.high_confidence_targets.startsWith('http://') || 
                            params.high_confidence_targets.startsWith('https://') || 
                            params.high_confidence_targets.startsWith('ftp://')
        ch_targets = Channel.from([
            ['high_confidence', file(params.high_confidence_targets, checkIfExists: !is_targets_remote)],
            ['gene_panel', file(params.gene_panel_targets, checkIfExists: !is_targets_remote)],
            ['wes_utr', file(params.wes_utr_targets, checkIfExists: !is_targets_remote)]
        ])
    } else {
        ch_targets = Channel.empty()
    }
    
    // Optional: ONT tandem repeats BED
    ch_tandem_repeats = params.tandem_repeats ? 
        Channel.value(file(params.tandem_repeats)) : 
        Channel.empty()
    
    emit:
    fasta            = ch_fasta           // channel: reference FASTA file
    fasta_fai        = ch_fasta_fai       // channel: reference FAI index
    benchmark_vcf    = ch_benchmark_vcf   // channel: truth VCF (may be empty)
    benchmark_vcf_tbi = ch_benchmark_vcf_tbi  // channel: truth VCF index (may be empty)
    targets          = ch_targets         // channel: [target_name, bed_file]
    tandem_repeats   = ch_tandem_repeats  // channel: tandem repeats BED (may be empty)
}
