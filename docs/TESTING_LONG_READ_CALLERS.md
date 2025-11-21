# Testing Long-Read SV Callers

This document describes how to test all structural variant callers in the pipeline, including PacBio and ONT long-read technologies.

## Overview

The pipeline supports multiple SV callers across all sequencing technologies:

### Illumina SV Callers
- **Manta**: Short-read structural variant caller

### PacBio HiFi SV Callers
- **CuteSV**: General-purpose SV caller optimized for long reads
- **PBSV**: PacBio-specific SV caller (two-stage: discover → call)

### Oxford Nanopore (ONT) SV Callers
- **CuteSV**: General-purpose SV caller optimized for long reads
- **Sniffles**: ONT-optimized SV caller with tandem repeat support

## Comprehensive Test Profile (`test_nfcore`)

The `test_nfcore` profile now tests **all sequencing technologies and SV callers** in a single pipeline run.

**Usage:**
```bash
nextflow run . -profile test_nfcore,docker
```

**What it tests:**
- ✅ Illumina WES BAM processing (Manta)
- ✅ Illumina WGS BAM processing (Manta)
- ✅ PacBio HiFi BAM processing (CuteSV + PBSV)
- ✅ Oxford Nanopore BAM processing (CuteSV + Sniffles)
- ✅ VCF compression and indexing for all callers
- ✅ Complete pipeline integration

**Configuration:**
- Resource limits: 2 CPUs, 6GB memory, 6h runtime
- Test data from nf-core/test-datasets:
  - Illumina: `test.paired_end.sorted.bam`, `test2.paired_end.sorted.bam`
  - PacBio: `test_hifi.sorted.bam`
  - ONT: `test.sorted.bam`
- Benchmarking is disabled (no SV truth set for test data)
- Increased time limit for long-read callers (2h per process)

## Running Tests Locally

### Prerequisites
- Nextflow (v23.04.0 or later)
- Docker or Singularity

### Quick Start

**Test all SV callers (all technologies):**
```bash
nextflow run . -profile test_nfcore,docker --outdir test_results
```

This single command tests:
- Illumina short-read callers (Manta)
- PacBio HiFi callers (CuteSV + PBSV)
- ONT callers (CuteSV + Sniffles)

## Expected Outputs

### Complete Test Outputs

```
test_results/
├── Illumina_WES/
│   └── Manta/
│       ├── Illumina_WES_Manta.vcf.gz
│       └── Illumina_WES_Manta.vcf.gz.tbi
├── Illumina_WGS/
│   └── Manta/
│       ├── Illumina_WGS_Manta.vcf.gz
│       └── Illumina_WGS_Manta.vcf.gz.tbi
├── PacBio/
│   ├── CuteSV/
│   │   ├── PacBio_CuteSV.vcf.gz
│   │   └── PacBio_CuteSV.vcf.gz.tbi
│   └── Pbsv/
│       ├── PacBio_Pbsv.vcf.gz
│       └── PacBio_Pbsv.vcf.gz.tbi
├── ONT/
│   ├── CuteSV/
│   │   ├── ONT_CuteSV.vcf.gz
│   │   └── ONT_CuteSV.vcf.gz.tbi
│   └── Sniffles/
│       ├── ONT_Sniffles.vcf.gz
│       └── ONT_Sniffles.vcf.gz.tbi
└── pipeline_info/
    ├── execution_timeline.html
    ├── execution_report.html
    ├── execution_trace.txt
    └── pipeline_dag.svg
```

**Total expected VCF files**: 6 callers × 2 files (vcf.gz + tbi) = 12 files

## CI/CD Testing

The GitHub Actions workflow automatically tests all SV callers on every push and pull request.

### CI Test Jobs

1. **profile-check**: Validates that test_nfcore profile loads correctly
2. **run-test**: Runs complete SV calling pipeline for all technologies
   - Tests Illumina (Manta)
   - Tests PacBio (CuteSV + PBSV)
   - Tests ONT (CuteSV + Sniffles)
   - Validates all VCF outputs

### Viewing CI Results

Check the Actions tab in the GitHub repository to see:
- Profile validation results
- Complete test execution logs
- Generated artifacts (all VCF files from all callers)
- Output validation summary showing which callers produced results

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
# Illumina
cat work/*/MANTA_GERMLINE/.command.log

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
nextflow run . -profile test_nfcore,docker -resume
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
- **Illumina WES**: `data/genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam`
- **Illumina WGS**: `data/genomics/homo_sapiens/illumina/bam/test2.paired_end.sorted.bam`
- **PacBio**: `data/genomics/homo_sapiens/pacbio/bam/test.sorted.bam`
- **ONT**: `data/genomics/homo_sapiens/nanopore/bam/test.sorted.bam`
- **Reference**: `data/genomics/homo_sapiens/genome/genome.fasta`

## Contributing

When adding new SV callers or test profiles:

1. Create or update test config file in `conf/test_nfcore.config`
2. Update the CI workflow (`.github/workflows/ci.yml`) if needed
3. Add expected outputs to output validation
4. Update this documentation with new caller information
5. Test locally with `nextflow run . -profile test_nfcore,docker`
6. Submit a pull request

## Related Documentation

- [Main Pipeline README](../README.md)
- [Configuration Profiles](../nextflow.config)
- [CI/CD Workflow](../.github/workflows/ci.yml)
