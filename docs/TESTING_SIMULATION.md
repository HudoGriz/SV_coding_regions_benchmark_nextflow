# Testing Guide: Simulation Features

This guide provides comprehensive instructions for testing the simulation and statistics features of the SV Coding Regions Benchmark pipeline.

## Table of Contents
- [Quick Start](#quick-start)
- [Test Configuration](#test-configuration)
- [Running Tests](#running-tests)
- [Validation](#validation)
- [CI/CD Integration](#cicd-integration)
- [Troubleshooting](#troubleshooting)
- [Advanced Testing](#advanced-testing)

## Quick Start

### Prerequisites
- Nextflow 24.04.0 or later
- Docker or Singularity
- Test data (included in `test_data/` directory)

### Run Full Simulation Test
```bash
# With Docker
nextflow run main.nf -profile test_simulation,docker

# With Singularity
nextflow run main.nf -profile test_simulation,singularity
```

**Expected Duration**: 5-10 minutes on standard hardware

## Test Configuration

### Test Profile: test_simulation

Located in `conf/test_simulation.config`, this profile includes:

```groovy
params {
    // Simulation enabled
    simulate_targets = true
    num_simulations = 5  // Small number for fast testing
    gencode_gtf = "${projectDir}/test_data/gencode_test.gtf.gz"
    
    // Statistics enabled
    gather_statistics = true
    
    // Single technology for speed
    illumina_wes_bam = "${projectDir}/test_data/illumina_wes.bam"
    
    // Resource limits
    max_cpus = 2
    max_memory = '6.GB'
    max_time = '6.h'
}
```

### Test Data

**gencode_test.gtf.gz**: Minimal GTF file with 3 genes
- TEST1 on chr1 (+ strand)
- TEST2 on chr1 (- strand)
- TEST3 on chr2 (+ strand)

Each gene has:
- Multiple exons
- CDS (coding sequence) regions
- UTR (untranslated regions)

## Running Tests

### 1. Full Pipeline Test

Tests all features: SV calling, benchmarking, simulation, and statistics.

```bash
nextflow run main.nf -profile test_simulation,docker
```

**What it does:**
1. Calls SVs with Manta on test Illumina WES data
2. Benchmarks against standard targets (high_conf, gene_panel, wes_utr)
3. Simulates 5 random target region sets
4. Benchmarks each simulated target set
5. Generates statistics and plots

### 2. Simulation Only Test

Test just the simulation module without SV calling.

```bash
nextflow run main.nf \
  -profile test_simulation,docker \
  --skip_benchmarking
```

### 3. Statistics Only Test

Test statistics generation on existing results.

```bash
nextflow run main.nf \
  -profile test_simulation,docker \
  --simulate_targets false \
  --gather_statistics true
```

### 4. Custom Number of Simulations

```bash
nextflow run main.nf \
  -profile test_simulation,docker \
  --num_simulations 10
```

### 5. With Custom GTF

```bash
nextflow run main.nf \
  -profile test_simulation,docker \
  --gencode_gtf /path/to/custom.gtf.gz \
  --num_simulations 20
```

### 6. Test with Reports

Generate detailed execution reports:

```bash
nextflow run main.nf \
  -profile test_simulation,docker \
  -with-report simulation_report.html \
  -with-timeline simulation_timeline.html \
  -with-trace simulation_trace.txt \
  -with-dag simulation_dag.html
```

## Validation

### Automated Validation Script

Create a validation script to check outputs:

```bash
#!/bin/bash
# validate_simulation.sh

OUTDIR="test_results_simulation"

echo "🔍 Validating simulation outputs..."

# Check simulated BED files
BED_COUNT=$(ls ${OUTDIR}/simulated_targets/simulation_*.bed 2>/dev/null | wc -l)
echo "✓ Found $BED_COUNT simulated BED files"

if [ "$BED_COUNT" -lt 5 ]; then
    echo "❌ Expected at least 5 BED files, found $BED_COUNT"
    exit 1
fi

# Check each BED file is valid
for bed in ${OUTDIR}/simulated_targets/simulation_*.bed; do
    if [ ! -s "$bed" ]; then
        echo "❌ Empty BED file: $bed"
        exit 1
    fi
    
    # Check BED format (at least 3 columns)
    COLS=$(head -1 "$bed" | awk '{print NF}')
    if [ "$COLS" -lt 3 ]; then
        echo "❌ Invalid BED format in $bed"
        exit 1
    fi
done
echo "✓ All BED files are valid"

# Check benchmarking results
if [ -d "${OUTDIR}/benchmarking" ]; then
    SUMMARY_COUNT=$(find ${OUTDIR}/benchmarking -name "summary.json" | wc -l)
    echo "✓ Found $SUMMARY_COUNT benchmark result summaries"
fi

# Check statistics outputs
if [ -d "${OUTDIR}/statistics" ]; then
    echo "✓ Statistics directory exists"
    
    # Check plots
    PLOT_COUNT=$(ls ${OUTDIR}/statistics/plots/*.png 2>/dev/null | wc -l)
    echo "✓ Found $PLOT_COUNT plot files"
    
    # Check tables
    TABLE_COUNT=$(ls ${OUTDIR}/statistics/tables/*.csv 2>/dev/null | wc -l)
    echo "✓ Found $TABLE_COUNT table files"
    
    # Check summary
    if [ -f "${OUTDIR}/statistics/summary_statistics.txt" ]; then
        echo "✓ Summary statistics file exists"
    fi
fi

echo "✅ Validation complete!"
```

Usage:
```bash
chmod +x validate_simulation.sh
./validate_simulation.sh
```

### Manual Validation Steps

#### 1. Verify Simulated BED Files

```bash
# Count files
ls test_results_simulation/simulated_targets/simulation_*.bed | wc -l

# Check first file
head test_results_simulation/simulated_targets/simulation_1.bed

# Verify BED format
awk '{print NF}' test_results_simulation/simulated_targets/simulation_1.bed | sort | uniq -c

# Check coverage
bedtools genomecov -i test_results_simulation/simulated_targets/simulation_1.bed \
  -g reference.fa.fai | head
```

#### 2. Verify Benchmarking Results

```bash
# List all benchmark results
find test_results_simulation/benchmarking -name "summary.json"

# Check a specific simulation result
cat test_results_simulation/benchmarking/illumina_wes/Manta/simulated/simulation_1/summary.json | jq .

# Compare metrics across simulations
for i in {1..5}; do
  echo "Simulation $i:"
  jq '.precision, .recall, .f1' \
    test_results_simulation/benchmarking/illumina_wes/Manta/simulated/simulation_${i}/summary.json
done
```

#### 3. Verify Statistics Outputs

```bash
# List generated plots
ls -lh test_results_simulation/statistics/plots/

# View summary statistics
cat test_results_simulation/statistics/summary_statistics.txt

# Check table format
head test_results_simulation/statistics/tables/*.csv
```

### Expected Output Structure

```
test_results_simulation/
├── manta/
│   └── illumina_wes/
│       └── illumina_wes.diploidSV.vcf.gz
├── benchmarking/
│   └── illumina_wes/
│       └── Manta/
│           ├── high_confidence/
│           │   └── summary.json
│           ├── gene_panel/
│           │   └── summary.json
│           ├── wes_utr/
│           │   └── summary.json
│           └── simulated/
│               ├── simulation_1/
│               │   └── summary.json
│               ├── simulation_2/
│               │   └── summary.json
│               ├── simulation_3/
│               │   └── summary.json
│               ├── simulation_4/
│               │   └── summary.json
│               └── simulation_5/
│                   └── summary.json
├── simulated_targets/
│   ├── simulation_1.bed
│   ├── simulation_2.bed
│   ├── simulation_3.bed
│   ├── simulation_4.bed
│   └── simulation_5.bed
└── statistics/
    ├── plots/
    │   ├── precision_recall.png
    │   ├── f1_scores.png
    │   ├── recall_by_caller.png
    │   └── precision_by_caller.png
    ├── tables/
    │   ├── summary_table.csv
    │   └── detailed_metrics.csv
    └── summary_statistics.txt
```

## CI/CD Integration

### GitHub Actions

The pipeline includes a comprehensive GitHub Actions workflow for testing simulation features.

**File**: `.github/workflows/test_simulation.yml`

**Triggers:**
- Push to master/main
- Pull requests to master/main
- Changes to simulation-related files
- Manual workflow dispatch

**Jobs:**
1. **test-simulation**: Full pipeline test with simulation
2. **test-simulation-modules**: Individual module testing
3. **test-workflow-syntax**: Syntax validation

**Running Manually:**
1. Go to GitHub Actions tab
2. Select "Test Simulation Features" workflow
3. Click "Run workflow"
4. Select branch and run

### Local CI Testing

Simulate CI environment locally using `act`:

```bash
# Install act
brew install act  # macOS
# or
curl https://raw.githubusercontent.com/nektos/act/master/install.sh | sudo bash

# Run simulation test workflow
act -j test-simulation

# Run with specific profile
act -j test-simulation -s GITHUB_TOKEN
```

## Troubleshooting

### Common Issues

#### 1. GTF File Not Found

**Error:**
```
ERROR ~ Error executing process > 'SIMULATE_TARGETS'
Caused by: Missing input file: gencode_test.gtf.gz
```

**Solution:**
```bash
# Check file exists
ls -lh test_data/gencode_test.gtf.gz

# If missing, recreate it
# (See test_data/README.md for instructions)
```

#### 2. Python Script Fails

**Error:**
```
ERROR ~ Error executing process > 'SIMULATE_TARGETS'
Python script failed with exit code 1
```

**Solution:**
```bash
# Check Python script syntax
python -m py_compile bin/Python/simulate_targets.py

# Test script manually
python bin/Python/simulate_targets.py --help

# Check Python version in container
docker run -it <container_image> python --version
```

#### 3. R Statistics Fail

**Error:**
```
ERROR ~ Error executing process > 'GATHER_STATISTICS'
R script failed
```

**Solution:**
```bash
# Check R scripts exist
ls -lh bin/R/*.R

# Verify R dependencies
docker run -it library://blazv/benchmark-sv/r-env:4-4-1 R

# In R, check packages:
# install.packages(c("ggplot2", "dplyr", "tidyr"))
```

#### 4. Insufficient Memory

**Error:**
```
ERROR ~ Process exceeded available memory
```

**Solution:**
```bash
# Increase memory in config
nextflow run main.nf \
  -profile test_simulation,docker \
  --max_memory 12.GB

# Or reduce simulation count
nextflow run main.nf \
  -profile test_simulation,docker \
  --num_simulations 3
```

#### 5. Empty Simulated BED Files

**Error:**
Simulation completes but BED files are empty

**Solution:**
```bash
# Check GTF has valid gene annotations
zcat test_data/gencode_test.gtf.gz | grep -c "gene_type"

# Verify GTF has exon/CDS features
zcat test_data/gencode_test.gtf.gz | cut -f3 | sort | uniq -c

# Check chromosome names match reference
zcat test_data/gencode_test.gtf.gz | cut -f1 | sort | uniq
```

### Debug Mode

Run with debug information:

```bash
nextflow run main.nf \
  -profile test_simulation,docker \
  -with-trace \
  -with-report \
  -with-timeline \
  -with-dag \
  --outdir test_debug \
  -ansi-log false
```

### Check Individual Processes

Test specific processes in isolation:

```bash
# Test simulation process only
nextflow run workflows/simulate_and_benchmark.nf \
  --fasta test_data/genome.fa \
  --gencode_gtf test_data/gencode_test.gtf.gz \
  --num_simulations 5

# Test statistics process only  
nextflow run workflows/analysis_and_plots.nf \
  --run_dir test_results_simulation
```

## Advanced Testing

### Stress Testing

Test with larger simulation counts:

```bash
# 100 simulations
nextflow run main.nf \
  -profile test_simulation,docker \
  --num_simulations 100 \
  --max_memory 16.GB \
  --max_cpus 8

# Measure performance
/usr/bin/time -v nextflow run main.nf -profile test_simulation,docker
```

### Parallel Execution

Test parallel simulation:

```bash
nextflow run main.nf \
  -profile test_simulation,docker \
  -process.executor local \
  -process.cpus 4 \
  --num_simulations 20
```

### Custom Test Data

Create custom test scenarios:

```bash
# Use your own GTF
nextflow run main.nf \
  -profile test_simulation,docker \
  --gencode_gtf /path/to/custom.gtf.gz \
  --num_simulations 10

# Use different genome regions
nextflow run main.nf \
  -profile test_simulation,docker \
  --high_confidence_targets custom_regions.bed
```

### Reproducibility Testing

Verify simulations are reproducible:

```bash
# Run 1
nextflow run main.nf -profile test_simulation,docker --outdir run1

# Run 2
nextflow run main.nf -profile test_simulation,docker --outdir run2

# Compare outputs
diff run1/simulated_targets/simulation_1.bed run2/simulated_targets/simulation_1.bed
```

### Performance Profiling

Profile resource usage:

```bash
nextflow run main.nf \
  -profile test_simulation,docker \
  -with-trace trace.txt

# Analyze resource usage
awk -F'\t' 'NR>1 {print $4, $5, $6}' trace.txt | \
  awk '{cpu+=$1; mem+=$2; time+=$3; n++} END {print "Avg CPU:", cpu/n, "Avg Mem:", mem/n, "Avg Time:", time/n}'
```

## Test Checklist

- [ ] Full simulation test passes
- [ ] All 5 BED files generated
- [ ] BED files have valid format
- [ ] Benchmarking completes for all simulations
- [ ] Statistics directory created
- [ ] Plots generated successfully
- [ ] Tables have expected format
- [ ] Summary statistics file created
- [ ] No error messages in log
- [ ] Resource usage within limits
- [ ] Output structure matches expected
- [ ] CI/CD workflow passes

## Resources

- [Nextflow Testing Best Practices](https://www.nextflow.io/docs/latest/testing.html)
- [nf-core Testing Guidelines](https://nf-co.re/docs/contributing/test_data)
- [GitHub Actions for Nextflow](https://github.com/nf-core/setup-nextflow)

## Support

For issues or questions:
1. Check this troubleshooting guide
2. Review test_data/README.md
3. Check GitHub Issues
4. Review Nextflow execution logs

## Contributing Tests

To add new test cases:

1. Create test data in `test_data/`
2. Add test configuration in `conf/`
3. Document in test README
4. Add CI workflow if needed
5. Update this testing guide

Example test contribution:
```bash
# 1. Create test data
cp my_test.gtf.gz test_data/

# 2. Create test config
cat > conf/test_my_feature.config << EOF
params {
    simulate_targets = true
    gencode_gtf = "${projectDir}/test_data/my_test.gtf.gz"
    # ... other params
}
EOF

# 3. Test it
nextflow run main.nf -profile test_my_feature,docker

# 4. Document and submit PR
```
