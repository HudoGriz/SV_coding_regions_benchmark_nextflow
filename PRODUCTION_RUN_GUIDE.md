# Production Run Guide

This guide explains how to run the Nextflow pipeline on your production data, replacing the original bash script approach.

## Overview

**Original approach:** `scripts/run_all_analyses_and_bench.sh`  
**New approach:** `nextflow run main.nf -params-file params_production.yaml`

## Quick Start

### 1. Update Paths in params_production.yaml

Edit `params_production.yaml` and update all file paths to match your actual data locations:

```yaml
# Required files:
fasta: '/path/to/reference/hs37d5.fa'
benchmark_vcf: '/path/to/HG002_SVs_Tier1_v0.6.vcf.gz'

# BAM files:
illumina_wes_bam: '/path/to/Illumina_wes.bam'
illumina_wgs_bam: '/path/to/Illumina_wgs.bam'
pacbio_bam: '/path/to/Pacbio.bam'
ont_bam: '/path/to/ONT.bam'

# Target regions:
high_confidence_targets: '/path/to/targets.bed'
# ... etc
```

### 2. Verify File Requirements

Make sure all required files exist and have proper indexes:

```bash
# BAM files need .bai indexes
ls -lh data/Illumina_wes/bam/*.bai
ls -lh data/Illumina_wgs/bam/*.bai
ls -lh data/Pacbio/bam/*.bai
ls -lh data/ONT/bam/*.bai

# Reference needs .fai index
ls -lh data/reference/*.fai

# Benchmark VCF needs .tbi index
ls -lh data/benchmark/*.tbi
```

If indexes are missing, create them:
```bash
# BAM indexes
samtools index file.bam

# FASTA index
samtools faidx reference.fa

# VCF index
tabix -p vcf file.vcf.gz
```

### 3. Choose Your Execution Method

#### Option A: Using the Helper Script (Recommended)

```bash
# 1. Make script executable
chmod +x run_production.sh

# 2. Edit paths in the script if needed
nano run_production.sh

# 3. Run
./run_production.sh
```

#### Option B: Direct Nextflow Command

```bash
nextflow run main.nf \
    -params-file params_production.yaml \
    -profile singularity \
    -resume
```

#### Option C: Override Specific Parameters

```bash
# Override individual parameters on the command line
nextflow run main.nf \
    -params-file params_production.yaml \
    -profile singularity \
    --skip_pbsv true \
    --max_cpus 16 \
    -resume
```

## Parameter Mapping

Here's how the old bash script parameters map to Nextflow:

| Bash Script | Nextflow Parameter |
|-------------|-------------------|
| `$project_dir` | Handled automatically by Nextflow |
| `run_name="benchmarking_run"` | `--run_name 'benchmarking_run'` |
| BAM file paths | `--illumina_wes_bam`, `--pacbio_bam`, etc. |
| `scripts/configs/config.sh` | Built into pipeline logic |
| Truvari params | `--truvari_refdist`, `--truvari_pctsize`, etc. |

## Output Structure

The pipeline produces the same output structure as the original bash script:

```
benchmarking_run/
├── calls/                      # SV calls (same as original)
│   ├── Illumina_wes/
│   │   └── sv/manta/
│   ├── Illumina_wgs/
│   │   └── sv/manta/
│   ├── Pacbio/
│   │   ├── sv/pbsv/
│   │   └── sv/cutesv/
│   └── ONT/
│       ├── sv/cutesv/
│       └── sv/sniffles/
├── real_intervals/             # Benchmarking results (same as original)
│   ├── Illumina_wes/
│   ├── Illumina_wgs/
│   ├── Pacbio/
│   └── ONT/
└── pipeline_info/              # NEW: Nextflow execution reports
    ├── report.html
    ├── timeline.html
    ├── trace.txt
    └── dag.html
```

## Key Differences from Bash Script

### Advantages of Nextflow:

1. **Resume capability**: Use `-resume` to continue from where it stopped
   ```bash
   # If pipeline fails, just run again with -resume
   nextflow run main.nf -params-file params_production.yaml -resume
   ```

2. **Parallel execution**: All technologies run in parallel automatically
   - Original: Sequential (Illumina WES → WGS → PacBio → ONT)
   - Nextflow: All run simultaneously (faster!)

3. **Better resource management**: Automatic CPU/memory allocation

4. **Detailed reports**: HTML reports with execution times, resource usage

5. **No manual tracking**: No need for `pipeline_status.log` file

6. **Containerization**: Each tool runs in its own container (reproducible!)

### Removed Manual Steps:

- ❌ No need for `is_done()` / `mark_done()` tracking
- ❌ No need for separate tool-specific scripts
- ❌ No need for manual config file management

## Troubleshooting

### Issue: "File not found" errors

**Solution:** Check that all paths in `params_production.yaml` are absolute paths

### Issue: "Permission denied" for containers

**Solution:** Make sure Singularity/Docker is properly configured:
```bash
# Test Singularity
singularity --version

# Test Docker
docker --version
```

### Issue: Memory errors

**Solution:** Adjust resource limits in `params_production.yaml`:
```yaml
max_cpus: 16      # Reduce if needed
max_memory: '64.GB'  # Reduce if needed
```

### Issue: Want to skip a technology

**Solution:** Don't provide the BAM file parameter:
```yaml
# Skip ONT by commenting out or removing:
# ont_bam: '/path/to/ont.bam'
```

### Issue: Want to skip PBSV (PacBio tool)

**Solution:** Set in params file:
```yaml
skip_pbsv: true
```

## Performance Comparison

### Original Bash Script (Sequential):
```
Illumina WES:  ~2 hours
Illumina WGS:  ~4 hours
PacBio:        ~3 hours
ONT:           ~3 hours
Benchmarking:  ~1 hour
----------------------------
Total:         ~13 hours
```

### Nextflow Pipeline (Parallel):
```
All technologies run simultaneously
Total:         ~4-5 hours (with sufficient resources)
```

## Advanced Usage

### Run only specific technologies

```bash
# Only Illumina WES
nextflow run main.nf \
    --illumina_wes_bam /path/to/wes.bam \
    --skip_benchmarking false \
    -profile singularity

# Only PacBio and ONT
nextflow run main.nf \
    --pacbio_bam /path/to/pacbio.bam \
    --ont_bam /path/to/ont.bam \
    --skip_benchmarking false \
    -profile singularity
```

### Use different Truvari parameters

```bash
nextflow run main.nf \
    -params-file params_production.yaml \
    --truvari_refdist 1000 \
    --truvari_pctsize 0.5 \
    -profile singularity
```

### Run with more/less resources

```bash
nextflow run main.nf \
    -params-file params_production.yaml \
    --max_cpus 32 \
    --max_memory '256.GB' \
    -profile singularity
```

## Need Help?

1. Check pipeline logs: `.nextflow.log`
2. Check work directory: `work/` contains all intermediate files
3. Check reports: `benchmarking_run/pipeline_info/`
4. Use `-with-tower` for cloud monitoring (if available)

## Validation

After running, verify outputs match original bash script:

```bash
# Check that all expected VCF files exist
ls -lh benchmarking_run/calls/*/sv/*/*.vcf.gz

# Check benchmarking results
ls -lh benchmarking_run/real_intervals/*/truvari/*/

# Compare file sizes (should be similar to original)
du -sh benchmarking_run/calls/
du -sh benchmarking_run/real_intervals/
```
