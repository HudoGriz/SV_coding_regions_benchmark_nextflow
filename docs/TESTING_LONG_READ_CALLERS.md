# Testing Long-Read SV Callers

This document describes how to test the PacBio and ONT structural variant callers in the pipeline.

## Overview

The pipeline supports multiple SV callers for long-read sequencing technologies:

### PacBio HiFi SV Callers
- **CuteSV**: General-purpose SV caller optimized for long reads
- **PBSV**: PacBio-specific SV caller (two-stage: discover → call)

### Oxford Nanopore (ONT) SV Callers
- **CuteSV**: General-purpose SV caller optimized for long reads
- **Sniffles**: ONT-optimized SV caller with tandem repeat support

## Test Profiles

### PacBio Test Profile (`test_pacbio`)

Tests both CuteSV and PBSV using nf-core test datasets.

**Usage:**
```bash
nextflow run . -profile test_pacbio,docker
```

**What it tests:**
- PacBio HiFi BAM processing
- CuteSV variant calling on PacBio data
- PBSV discover and call workflow
- VCF compression and indexing
- Output file generation

**Configuration:**
- Resource limits: 2 CPUs, 6GB memory, 6h runtime
- Uses `test_hifi.sorted.bam` from nf-core test datasets
- Benchmarking is disabled (no SV truth set for test data)
- Error strategy: ignore (to identify which callers work)

### ONT Test Profile (`test_ont`)

Tests both CuteSV and Sniffles using nf-core test datasets.

**Usage:**
```bash
nextflow run . -profile test_ont,docker
```

**What it tests:**
- Oxford Nanopore BAM processing
- CuteSV variant calling on ONT data
- Sniffles variant calling with tandem repeat annotation
- VCF compression and indexing
- Output file generation

**Configuration:**
- Resource limits: 2 CPUs, 6GB memory, 6h runtime
- Uses `test.sorted.bam` from nf-core nanopore test datasets
- Benchmarking is disabled (no SV truth set for test data)
- Error strategy: ignore (to identify which callers work)

## Running Tests Locally

### Prerequisites
- Nextflow (v23.04.0 or later)
- Docker or Singularity

### Quick Start

**Test PacBio callers:**
```bash
nextflow run . -profile test_pacbio,docker --outdir results_pacbio
```

**Test ONT callers:**
```bash
nextflow run . -profile test_ont,docker --outdir results_ont
```

**Test all long-read callers:**
```bash
# PacBio
nextflow run . -profile test_pacbio,docker --outdir results_pacbio

# ONT
nextflow run . -profile test_ont,docker --outdir results_ont
```

## Expected Outputs

### PacBio Test Outputs

```
results_pacbio/
├── PacBio/
│   ├── CuteSV/
│   │   ├── PacBio_CuteSV.vcf.gz
│   │   └── PacBio_CuteSV.vcf.gz.tbi
│   └── Pbsv/
│       ├── PacBio_Pbsv.vcf.gz
│       └── PacBio_Pbsv.vcf.gz.tbi
└── pipeline_info/
```

### ONT Test Outputs

```
results_ont/
├── ONT/
│   ├── CuteSV/
│   │   ├── ONT_CuteSV.vcf.gz
│   │   └── ONT_CuteSV.vcf.gz.tbi
│   └── Sniffles/
│       ├── ONT_Sniffles.vcf.gz
│       └── ONT_Sniffles.vcf.gz.tbi
└── pipeline_info/
```

## CI/CD Testing

The GitHub Actions workflow automatically tests both profiles on every push and pull request.

### CI Test Jobs

1. **profile-check**: Validates that test_pacbio and test_ont profiles load correctly
2. **test-pacbio**: Runs full PacBio SV calling pipeline
3. **test-ont**: Runs full ONT SV calling pipeline

### Viewing CI Results

Check the Actions tab in the GitHub repository to see:
- Profile validation results
- Test execution logs
- Generated artifacts (VCF files, logs)

## Troubleshooting

### Common Issues

**Issue**: Test data download fails
- **Solution**: Check internet connectivity and nf-core test-datasets availability

**Issue**: Docker/Singularity not found
- **Solution**: Install container engine or use a different profile

**Issue**: Out of memory errors
- **Solution**: Increase memory limits in the test config files

**Issue**: Caller produces no output
- **Solution**: Check `.nextflow.log` for detailed error messages. The error strategy is set to 'ignore' to allow partial completion.

### Debugging

**View detailed logs:**
```bash
cat .nextflow.log
```

**Check specific process logs:**
```bash
# PacBio
cat work/*/CUTESV_PACBIO/.command.log
cat work/*/PBSV_DISCOVER/.command.log
cat work/*/PBSV_CALL/.command.log

# ONT
cat work/*/CUTESV_ONT/.command.log
cat work/*/SNIFFLES/.command.log
```

**Resume failed runs:**
```bash
nextflow run . -profile test_pacbio,docker -resume
```

## Customizing Tests

### Using Custom Test Data

Create a custom config file:

```groovy
// custom_test.config
params {
    pacbio_bam = '/path/to/your/pacbio.bam'
    ont_bam = '/path/to/your/ont.bam'
    fasta = '/path/to/reference.fasta'
    
    // Use test regions
    high_confidence_targets = "${test_data_base}/genomics/homo_sapiens/genome/genome.bed"
    gene_panel_targets = "${test_data_base}/genomics/homo_sapiens/genome/genome.bed"
    wes_utr_targets = "${test_data_base}/genomics/homo_sapiens/genome/genome.bed"
    
    outdir = 'custom_test_results'
}
```

Run with:
```bash
nextflow run . -profile docker -c custom_test.config
```

### Modifying Caller Parameters

Edit `conf/modules.config` to customize SV caller parameters:

```groovy
process {
    withName: 'CUTESV_PACBIO' {
        ext.args = '--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9'
    }
    
    withName: 'SNIFFLES' {
        ext.args = '--minsvlen 50 --minsupport 2'
    }
}
```

## Test Data Sources

All test data comes from the nf-core test-datasets repository:
- **Repository**: https://github.com/nf-core/test-datasets
- **Branch**: modules
- **PacBio data**: `data/genomics/homo_sapiens/pacbio/bam/test_hifi.sorted.bam`
- **ONT data**: `data/genomics/homo_sapiens/nanopore/bam/test.sorted.bam`
- **Reference**: `data/genomics/homo_sapiens/genome/genome.fasta`

## Contributing

When adding new callers or test profiles:

1. Create a new test config file in `conf/`
2. Add the profile to the CI workflow matrix
3. Create a dedicated test job if needed
4. Update this documentation
5. Submit a pull request

## Related Documentation

- [Main Pipeline README](../README.md)
- [Configuration Profiles](../nextflow.config)
- [CI/CD Workflow](../.github/workflows/ci.yml)
