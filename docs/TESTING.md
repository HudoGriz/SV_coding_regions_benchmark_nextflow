# Testing Guide for SV Calling Pipeline

This document describes how to test the SV calling pipeline both locally and via GitHub Actions CI/CD.

## Overview

The pipeline includes automated tests that:
- ✅ Validate Nextflow DSL2 syntax
- ✅ Test workflow logic with synthetic data
- ✅ Run in stub mode to verify process connections
- ✅ Check configuration validity
- ✅ Ensure compatibility with multiple Nextflow versions

## Test Data

### Synthetic Test Data

The tests use minimal synthetic data to keep CI fast:

| Component | Description | Size |
|-----------|-------------|------|
| Reference genome | Synthetic chr22 (1 MB) | ~1 MB |
| BAM files | 20 aligned reads | ~1 KB |
| Truth VCF | 5 synthetic SVs (DEL, INS, DUP) | ~1 KB |
| Target BED files | 3 region sets | ~1 KB |

**Total test data size:** ~1 MB (downloads quickly in CI)

## Local Testing

### Prerequisites

- Nextflow >= 23.04.0
- samtools
- bgzip and tabix (from htslib)
- Python 3
- Docker or Singularity (for full run, not required for stub mode)

### Quick Test (Stub Mode)

The fastest way to test the pipeline logic without executing actual processes:

```bash
# Run the automated test script
./test_local.sh
```

This script:
1. Creates synthetic test data
2. Generates test parameters
3. Runs the pipeline in stub mode
4. Validates workflow logic

**Stub mode** creates empty/mock outputs but verifies:
- Process connections are correct
- Channel operations work
- No syntax errors
- Configuration is valid

### Full Local Test (With Containers)

To run the full pipeline with actual tool execution:

```bash
# Generate test data (if not already done)
./test_local.sh

# Run with Docker
nextflow run main.nf -params-file test_params.yaml -profile test

# OR run with Singularity
nextflow run main.nf -params-file test_params.yaml -profile test,singularity
```

**Note:** Full execution requires downloading Docker/Singularity containers (~1-5 GB depending on tools).

### Manual Test Data Creation

If you want to create test data manually:

```bash
# Create directories
mkdir -p test_data

# Create synthetic reference (1 MB, chr22)
echo ">chr22" > test_data/reference.fa
python3 -c "import random; random.seed(42); print(''.join(random.choices('ACGT', k=1000000)))" >> test_data/reference.fa
samtools faidx test_data/reference.fa

# Create minimal BAM
cat > test_data/test.sam << EOF
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr22	LN:1000000
@RG	ID:test	SM:HG002	PL:ILLUMINA
read_1	0	chr22	10000	60	100M	*	0	0	ACGT...	*
EOF

samtools view -Sb test_data/test.sam | samtools sort -o test_data/test_sorted.bam
samtools index test_data/test_sorted.bam

# Create truth VCF with synthetic SVs
cat > test_data/truth.vcf << EOF
##fileformat=VCFv4.2
##contig=<ID=chr22,length=1000000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG002
chr22	10000	sv1	N	<DEL>	60	PASS	SVTYPE=DEL;SVLEN=-500	GT	0/1
EOF

bgzip -c test_data/truth.vcf > test_data/truth.vcf.gz
tabix -p vcf test_data/truth.vcf.gz

# Create target BED files
echo -e "chr22\t0\t1000000" > test_data/high_confidence.bed
echo -e "chr22\t10000\t100000" > test_data/gene_panel.bed
echo -e "chr22\t5000\t15000" > test_data/wes_utr.bed
```

### Cleaning Up

Remove all test files:

```bash
rm -rf test_data/ test_results/ test_params.yaml work/ .nextflow*
```

## GitHub Actions CI

### Automated Testing

The CI workflow (`.github/workflows/ci.yml`) runs automatically on:
- Push to `master`, `main`, or `dev` branches
- Pull requests to these branches
- Manual trigger via GitHub Actions UI

### CI Workflow Jobs

#### 1. Test Job

**Duration:** ~3-5 minutes

**Steps:**
1. Checkout code
2. Install Nextflow (latest stable)
3. Set up Apptainer/Singularity
4. Generate synthetic test data
5. Create test configuration
6. Run pipeline in stub mode
7. Validate outputs
8. Upload test artifacts

**What it validates:**
- ✅ Nextflow syntax is correct
- ✅ Workflow logic is sound
- ✅ Process connections work
- ✅ Configuration is valid
- ✅ Pipeline can be instantiated

#### 2. Lint Job

**Duration:** ~1-2 minutes

**Steps:**
1. Checkout code
2. Install Nextflow
3. Validate configuration syntax
4. Check DSL2 syntax
5. Display available profiles

**What it validates:**
- ✅ Config file syntax
- ✅ Profile definitions
- ✅ No obvious syntax errors

### Viewing CI Results

1. Go to your repository on GitHub
2. Click the "Actions" tab
3. Select the workflow run
4. View job logs and artifacts

### CI Artifacts

Each test run uploads:
- Test results directory
- Nextflow logs (`.nextflow.log`)
- Retention: 7 days

### Triggering Manual Runs

1. Go to Actions tab
2. Select "CI Tests" workflow
3. Click "Run workflow"
4. Select branch
5. Click "Run workflow" button

## Test Profile

The `test` profile (in `nextflow.config`) is optimized for CI/CD:

```groovy
test {
    // Minimal resources for CI
    params.max_cpus = 2
    params.max_memory = '6.GB'
    params.max_time = '10.m'
    
    // Use Docker for GitHub Actions
    docker.enabled = true
    
    // Override process resources
    process {
        cpus = 1
        memory = 2.GB
        time = 5.m
    }
}
```

**Key features:**
- Reduced resource requirements
- Docker-based (better for CI than Singularity)
- Lenient caching
- Fast timeouts to catch hanging processes

## Troubleshooting

### Test Failures

#### Syntax Errors

```
ERROR ~ No such variable: someVar
```

**Fix:** Check variable names in your Nextflow scripts. Ensure all variables are properly defined.

#### Missing Dependencies

```
ERROR ~ Cannot find samtools
```

**Fix:** Ensure all required tools are in PATH or update the CI workflow to install them.

#### Container Issues

```
ERROR ~ Unable to pull Docker image
```

**Fix:** 
- Check Docker Hub availability
- Verify container names in config
- Consider using alternative registries

#### Memory/Timeout Errors

```
ERROR ~ Process exceeded memory limit
```

**Fix:**
- Reduce test data size
- Increase test profile resources
- Check for memory leaks in processes

### Local Test Failures

#### Stub Mode Works But Full Run Fails

This usually indicates:
- Container issues (missing/incorrect containers)
- Tool-specific errors
- Resource constraints

**Debug:**
```bash
# Run with verbose logging
nextflow run main.nf -params-file test_params.yaml -profile test -with-trace -with-report

# Check process logs
cat work/*/*/.command.log
```

#### BAM Index Issues

```
ERROR ~ BAM index not found
```

**Fix:**
```bash
samtools index test_data/test_sorted.bam
```

## Advanced Testing

### Testing Specific Components

Test only Illumina WES:

```yaml
# test_params_wes.yaml
illumina_wes_bam: 'test_data/test_sorted.bam'
# Comment out other BAM parameters
```

```bash
nextflow run main.nf -params-file test_params_wes.yaml -profile test -stub-run
```

### Testing with Real (Small) Data

For more realistic testing, use small chromosomes or regions:

```bash
# Download small region (e.g., chr22:10000000-11000000)
samtools view -b ftp://path/to/full.bam chr22:10000000-11000000 > test_data/chr22_1MB.bam
samtools index test_data/chr22_1MB.bam

# Extract corresponding reference region
samtools faidx reference.fa chr22:10000000-11000000 > test_data/chr22_1MB.fa
```

### Performance Testing

Measure pipeline performance:

```bash
nextflow run main.nf \
  -params-file test_params.yaml \
  -profile test \
  -with-trace \
  -with-timeline \
  -with-report \
  -with-dag
```

Outputs:
- `trace.txt` - Process execution details
- `timeline.html` - Execution timeline
- `report.html` - Resource usage report
- `dag.svg` - Pipeline DAG visualization

## Best Practices

1. **Run stub mode first** - Catches logic errors quickly
2. **Test locally before pushing** - Saves CI time
3. **Use small data** - Keeps tests fast
4. **Test on clean environment** - Catches missing dependencies
5. **Check all profiles** - Ensure configs work across environments
6. **Document test cases** - Help future contributors

## CI/CD Badge

Add to your README.md:

```markdown
[![CI Tests](https://github.com/YOUR_USERNAME/YOUR_REPO/actions/workflows/ci.yml/badge.svg)](https://github.com/YOUR_USERNAME/YOUR_REPO/actions/workflows/ci.yml)
```

Replace `YOUR_USERNAME` and `YOUR_REPO` with your actual repository details.

## Future Enhancements

Potential test improvements:

- [ ] Add integration tests with real (small) BAM files
- [ ] Test multiple genome builds (GRCh37/38)
- [ ] Add parameter validation tests
- [ ] Test resume functionality
- [ ] Add performance regression tests
- [ ] Test cloud execution (AWS Batch, Google Batch)
- [ ] Add output validation checks
- [ ] Test error handling and edge cases

## Support

If you encounter testing issues:

1. Check the [Nextflow documentation](https://www.nextflow.io/docs/latest/)
2. Review CI logs carefully
3. Test locally with `-with-trace` and `-with-report`
4. Open an issue with:
   - Nextflow version
   - Error message
   - Relevant logs
   - Steps to reproduce

## Summary

Testing your pipeline is essential for:
- 🐛 Catching bugs early
- 🔄 Ensuring reproducibility
- 📦 Validating refactoring
- 🚀 Confident deployment
- 👥 Easier collaboration

The provided test setup gives you a solid foundation for continuous integration and quality assurance!
