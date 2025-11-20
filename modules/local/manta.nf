process MANTA {
    tag "${meta.id}_${meta.tool}"
    label 'process_medium'
    
    publishDir "${params.outdir}/calls/${meta.technology}/sv/manta", mode: 'copy'
    
    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_fai
    
    output:
    tuple val(meta), path("results/variants/diploidSV.vcf.gz"), path("results/variants/diploidSV.vcf.gz.tbi"), emit: vcf
    path "results/**", emit: all_results
    
    script:
    """
    # Step 1: Configure Manta
    configManta.py \\
        --bam ${bam} \\
        --referenceFasta ${fasta} \\
        --runDir manta_run
    
    # Step 2: Run Manta workflow
    manta_run/runWorkflow.py -j ${task.cpus}
    
    # Rename output directory for publishing
    mv manta_run results
    """
}
