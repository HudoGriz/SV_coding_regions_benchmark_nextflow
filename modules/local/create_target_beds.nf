process CREATE_EXOME_UTR_BED {
    tag "${genome_build}"
    label 'process_low'
    
    publishDir "${params.project_dir}/data/references", mode: 'copy'
    
    input:
    path gtf
    path r_container
    val genome_build
    
    output:
    path "exome_utr_*.bed", emit: bed
    
    script:
    def output_name = genome_build == 'GRCh38' ? 
        'exome_utr_gtf_GRCh38.bed' : 'exome_utr_gtf.bed'
    def r_script = genome_build == 'GRCh38' ? 
        'create_gencode_target_bed_GRCh38.R' : 'create_gencode_target_bed.R'
    """
    # Run R script to create exome + UTR BED file
    singularity exec \\
        -B \$(pwd) \\
        -B ${params.project_dir} \\
        ${r_container} \\
        Rscript ${params.project_dir}/scripts/R/${r_script}
    
    # Rename output if needed
    if [ -f "exome_utr_gtf.bed" ] && [ "${output_name}" != "exome_utr_gtf.bed" ]; then
        mv exome_utr_gtf.bed ${output_name}
    fi
    """
    
    stub:
    def output_name = genome_build == 'GRCh38' ? 
        'exome_utr_gtf_GRCh38.bed' : 'exome_utr_gtf.bed'
    """
    touch ${output_name}
    """
}
