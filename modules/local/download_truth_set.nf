process DOWNLOAD_GIAB_TRUTH_SET {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${output_dir}", mode: 'copy'
    
    input:
    tuple val(meta), val(base_url), val(output_dir)
    
    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "*.bed", emit: bed
    
    script:
    """
    # Download GIAB SV truth set files
    wget ${base_url}/HG002_SVs_Tier1_v0.6.bed
    wget ${base_url}/HG002_SVs_Tier1_v0.6.vcf.gz
    
    # Index VCF file
    tabix -p vcf HG002_SVs_Tier1_v0.6.vcf.gz
    """
    
    stub:
    """
    touch HG002_SVs_Tier1_v0.6.bed
    touch HG002_SVs_Tier1_v0.6.vcf.gz
    touch HG002_SVs_Tier1_v0.6.vcf.gz.tbi
    """
}

process DOWNLOAD_GIAB_TRUTH_SET_GRCH38 {
    tag "${meta.id}"
    label 'process_low'
    
    publishDir "${output_dir}", mode: 'copy'
    
    input:
    tuple val(meta), val(base_url), val(output_dir), val(genome)
    
    output:
    tuple val(meta), path("*_stvar.vcf.gz"), path("*_stvar.vcf.gz.tbi"), emit: vcf
    path "*_stvar.benchmark.bed", emit: bed
    
    script:
    """
    # Download GIAB SV truth set files for GRCh38
    wget ${base_url}/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed
    wget ${base_url}/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz
    wget ${base_url}/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz.tbi
    """
    
    stub:
    """
    touch GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed
    touch GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz
    touch GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz.tbi
    """
}

process DOWNLOAD_GIAB_TRUTH_SET_GRCH37_LIFTOVER {
    tag "${meta.id}"
    label 'process_low'
    
    publishDir "${output_dir}", mode: 'copy'
    
    input:
    tuple val(meta), val(base_url), val(output_dir), val(genome)
    
    output:
    tuple val(meta), path("*_stvar.vcf.gz"), path("*_stvar.vcf.gz.tbi"), emit: vcf
    path "*_stvar.benchmark.bed", emit: bed
    
    script:
    """
    # Download GIAB SV truth set files for GRCh37 (liftover from GRCh38)
    wget ${base_url}/GRCh37_HG002-T2TQ100-V1.0_stvar.benchmark.bed
    wget ${base_url}/GRCh37_HG002-T2TQ100-V1.0_stvar.vcf.gz
    wget ${base_url}/GRCh37_HG002-T2TQ100-V1.0_stvar.vcf.gz.tbi
    """
    
    stub:
    """
    touch GRCh37_HG002-T2TQ100-V1.0_stvar.benchmark.bed
    touch GRCh37_HG002-T2TQ100-V1.0_stvar.vcf.gz
    touch GRCh37_HG002-T2TQ100-V1.0_stvar.vcf.gz.tbi
    """
}
