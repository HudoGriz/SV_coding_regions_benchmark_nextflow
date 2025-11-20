process DOWNLOAD_SINGULARITY_IMAGES {
    tag "singularity_images"
    label 'process_low'
    
    publishDir "${params.project_dir}/singularity_images", mode: 'copy'
    
    input:
    val output_dir
    
    output:
    path "*.sif", emit: images
    path "r-env_4-4-1.sif", emit: r_env
    
    script:
    """
    # Download all required Singularity images
    singularity pull manta_latest.sif docker://dceoy/manta:latest
    singularity pull samtools_latest.sif docker://quay.io/biocontainers/samtools:1.19--h50ea8bc_0
    singularity pull cutesv_latest.sif docker://quay.io/biocontainers/cutesv:2.1.1--pyhdfd78af_0
    singularity pull pbsv_latest.sif docker://quay.io/pacbio/pbsv:2.10.0_build1
    singularity pull sniffles_latest.sif docker://quay.io/biocontainers/sniffles:2.4--pyhdfd78af_0
    singularity pull bedtools_latest.sif docker://quay.io/biocontainers/bedtools:2.31.1--h13024bc_3
    singularity pull truvari_modded.sif library://blazv/benchmark-sv/truvari_modded:latest
    singularity pull r-env_4-4-1.sif library://blazv/benchmark-sv/r-env:4-4-1
    """
    
    stub:
    """
    touch manta_latest.sif
    touch samtools_latest.sif
    touch cutesv_latest.sif
    touch pbsv_latest.sif
    touch sniffles_latest.sif
    touch bedtools_latest.sif
    touch truvari_modded.sif
    touch r-env_4-4-1.sif
    """
}
