# Guide: Preparing GIAB Resources

This guide explains how to use the GIAB resource preparation workflow to download and prepare validation resources for your pipeline.

## Overview

The `PREPARE_GIAB_RESOURCES` workflow automates the download and preparation of:
- GIAB HG002 structural variant truth sets
- Tandem repeat annotations
- GENCODE gene annotations
- Exome + UTR target BED files

## Quick Start

### Using Profiles (Recommended)

The easiest way to prepare GIAB resources is using the built-in profiles:

**For GRCh37/hs37d5:**
```bash
nextflow run main.nf -profile giab_grch37
```

**For GRCh38:**
```bash
nextflow run main.nf -profile giab_grch38
```

### Manual Configuration

You can also enable resource preparation manually:

```bash
nextflow run main.nf \
    --prepare_giab_resources \
    --genome hs37d5 \
    --project_dir /path/to/project \
    --r_container /path/to/r_container.sif
```

## Output Structure

### GRCh37/hs37d5

```
project_dir/
├── data/
│   ├── HG002_references/
│   │   ├── HG002_SVs_Tier1_v0.6.vcf.gz
│   │   ├── HG002_SVs_Tier1_v0.6.vcf.gz.tbi
│   │   ├── HG002_SVs_Tier1_v0.6.bed
│   │   └── human_hs37d5.trf.bed
│   └── references/
│       ├── gencode.v19.annotation.gtf.gz
│       └── exome_utr_gtf.bed
```

### GRCh38

```
project_dir/
├── data/
│   ├── HG002_references/
│   │   ├── GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz
│   │   ├── GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz.tbi
│   │   ├── GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed
│   │   └── human_GRCh38_no_alt_analysis_set.trf.bed
│   └── references/
│       ├── gencode.v49.annotation.gtf.gz
│       └── exome_utr_gtf_GRCh38.bed
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `prepare_giab_resources` | `false` | Enable GIAB resource preparation |
| `genome` | `'hs37d5'` | Genome build (`'hs37d5'` or `'GRCh38'`) |
| `project_dir` | `${projectDir}` | Project root directory |
| `r_container` | `${projectDir}/singularity_images/r_container.sif` | Path to R Singularity container |

## Prerequisites

### Required Software

- `wget` - for downloading files
- `tabix` - for VCF indexing (GRCh37 only)
- Singularity - for containerized execution

### Required Containers

You need an R container with the following packages:
- `rtracklayer`
- `GenomicFeatures`
- `dplyr`

Example container definition:
```singularity
Bootstrap: docker
From: rocker/tidyverse:latest

%post
    R -e 'install.packages("BiocManager")'
    R -e 'BiocManager::install(c("rtracklayer", "GenomicFeatures"))'
```

## Using Resources in Your Pipeline

After preparing resources, update your pipeline configuration:

**For GRCh37:**
```bash
nextflow run main.nf \
    --fasta /path/to/hs37d5.fa \
    --benchmark_vcf ${project_dir}/data/HG002_references/HG002_SVs_Tier1_v0.6.vcf.gz \
    --high_confidence_targets ${project_dir}/data/HG002_references/HG002_SVs_Tier1_v0.6.bed \
    --wes_utr_targets ${project_dir}/data/references/exome_utr_gtf.bed \
    --tandem_repeats ${project_dir}/data/HG002_references/human_hs37d5.trf.bed \
    # ... other parameters
```

**For GRCh38:**
```bash
nextflow run main.nf \
    --fasta /path/to/GRCh38.fa \
    --benchmark_vcf ${project_dir}/data/HG002_references/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz \
    --high_confidence_targets ${project_dir}/data/HG002_references/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed \
    --wes_utr_targets ${project_dir}/data/references/exome_utr_gtf_GRCh38.bed \
    --tandem_repeats ${project_dir}/data/HG002_references/human_GRCh38_no_alt_analysis_set.trf.bed \
    # ... other parameters
```

## Data Sources

### GIAB Truth Sets

**GRCh37/hs37d5:**
- Source: NCBI FTP GIAB Release
- Version: v3.3.2
- URL: ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37

**GRCh38:**
- Source: NCBI FTP GIAB Release
- Version: CMRG v1.00 (T2TQ100-V1.0)
- URL: ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38

### Tandem Repeats

- Source: speedseq annotations (Hall Lab)
- Repository: github.com/hall-lab/speedseq

### GENCODE Annotations

**GRCh37:**
- Version: v19 (Ensembl 75)
- URL: ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19

**GRCh38:**
- Version: v49
- URL: ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49

## Troubleshooting

### Download Failures

If downloads fail, the workflow will automatically retry up to 3 times. Check:
- Internet connectivity
- FTP server availability
- Firewall settings

### R Script Errors

If the exome+UTR BED creation fails:
- Verify R container has all required packages
- Check GTF file downloaded correctly
- Ensure sufficient memory (8GB recommended)

### File Permission Errors

Ensure the project directory is writable:
```bash
chmod -R u+w ${project_dir}/data
```

## Advanced Usage

### Custom URLs

You can modify the URLs in the workflow file if needed:

```groovy
// In workflows/prepare_giab_resources.nf
def giab_base_url = 'https://your-mirror.com/giab/data'
```

### Skip Specific Steps

To download only specific resources, you can comment out sections in the workflow or create a custom version.

### Different GIAB Samples

Currently configured for HG002. To use different samples (HG003, HG004, etc.), modify:
- Sample ID in metadata
- URLs to appropriate GIAB releases

## Integration with Existing Workflows

The resource preparation workflow is independent and can be run separately:

```bash
# Step 1: Prepare resources
nextflow run main.nf -profile giab_grch37

# Step 2: Run main pipeline (later)
nextflow run main.nf \
    --illumina_wes_bam sample.bam \
    --benchmark_vcf data/HG002_references/HG002_SVs_Tier1_v0.6.vcf.gz \
    # ... other parameters
```

## Notes

- Resource preparation is idempotent - safe to run multiple times
- Downloaded files are stored in publishDir, not work directory
- Total download size: ~500MB for GRCh37, ~800MB for GRCh38
- Expected runtime: 10-30 minutes depending on network speed
