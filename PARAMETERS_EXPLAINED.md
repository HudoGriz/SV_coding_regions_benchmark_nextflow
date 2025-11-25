# Complete Parameters Guide

This document explains **every parameter** you can use in the pipeline.

## 🔴 REQUIRED Parameters (Must Have)

### Absolutely Required
```yaml
fasta: '/path/to/reference.fa'
# Reference genome FASTA file (hs37d5, GRCh37, or GRCh38)
# Index (.fai) will be created automatically if missing

high_confidence_targets: '/path/to/targets.bed'
gene_panel_targets: '/path/to/gene_panel.bed'
wes_utr_targets: '/path/to/wes_utr.bed'
# Target BED files - ALWAYS required even if skip_benchmarking=true
# Can be dummy files (single line BED) if not doing benchmarking
```

### At Least One Required
```yaml
illumina_wes_bam: '/path/to/illumina_wes.bam'  # Illumina WES
illumina_wgs_bam: '/path/to/illumina_wgs.bam'  # Illumina WGS
pacbio_bam: '/path/to/pacbio.bam'              # PacBio HiFi
ont_bam: '/path/to/ont.bam'                    # Oxford Nanopore
# Provide at least ONE BAM file
# Each BAM must have a .bai index in the same directory
```

---

## 🟡 REQUIRED for Benchmarking

If `skip_benchmarking: false`, you must also provide:

```yaml
benchmark_vcf: '/path/to/truth.vcf.gz'
# Truth set VCF (e.g., GIAB HG002 Tier 1 v0.6)
# Must be bgzipped and have .tbi index
```

---

## 🟢 OPTIONAL Parameters (Have Defaults)

### Run Configuration
```yaml
run_name: 'benchmarking_run'
# Name for this run (default: 'benchmarking_run')

outdir: '/path/to/output'
# Output directory (default: ${launchDir}/${run_name})
```

### Pipeline Control
```yaml
skip_benchmarking: false
# Skip Truvari benchmarking (default: false)
# Set to true if you only want SV calls, not benchmarking

skip_pbsv: false
# Skip PBSV for PacBio data (default: false)
# If true, only CuteSV will run for PacBio

prepare_giab_resources: false
# Download and prepare GIAB resources (default: false)
# Special mode to download benchmark files

prepare_complete_data: false
# Download complete dataset including BAMs (default: false)
# Special mode to download all data
```

### Truvari Benchmarking Parameters
```yaml
truvari_refdist: 500
# Reference distance for matching variants (default: 500)

truvari_pctsize: 0.7
# Percent size similarity threshold (default: 0.7)

truvari_pctovl: 0.0
# Percent overlap threshold (default: 0.0)

truvari_pctseq: 0.0
# Percent sequence similarity threshold (default: 0.0)
```

### Simulation and Analysis
```yaml
simulate_targets: false
# Enable target region simulation (default: false)

num_simulations: 100
# Number of simulations to run (default: 100)

gencode_gtf: null
# GENCODE GTF annotation file for simulation
# Required if simulate_targets=true

gather_statistics: false
# Generate statistics and plots (default: false)
```

### Optional Target Files
```yaml
tandem_repeats: null
# Tandem repeats BED file for Sniffles annotation
# Optional, only used if running ONT data
```

### Resource Limits
```yaml
max_cpus: 24
# Maximum CPUs to use per process (default: 24)

max_memory: '128.GB'
# Maximum memory to allocate (default: 128.GB)

max_time: '48.h'
# Maximum time per process (default: 48.h)
```

### Advanced Data Preparation
```yaml
genome: 'hs37d5'
# Genome version for data preparation (default: 'hs37d5')
# Options: 'hs37d5', 'GRCh37', 'GRCh38'

skip_singularity_download: false
# Skip downloading Singularity containers (default: false)

skip_bam_download: false
# Skip downloading BAM files (default: false)

skip_reference_download: false
# Skip downloading reference genome (default: false)

download_grch37_liftover: false
# Download GRCh37 liftover files (default: false)
# Only for GRCh38 preparation mode
```

---

## 📊 Parameter Groups by Use Case

### Use Case 1: Minimal SV Calling (No Benchmarking)

```yaml
# REQUIRED
fasta: '/path/to/reference.fa'
illumina_wes_bam: '/path/to/sample.bam'
high_confidence_targets: '/path/to/dummy.bed'
gene_panel_targets: '/path/to/dummy.bed'
wes_utr_targets: '/path/to/dummy.bed'

# RECOMMENDED
skip_benchmarking: true
run_name: 'my_sv_calling'
max_cpus: 8
max_memory: '32.GB'
```

**Purpose:** Just call SVs, no comparison to truth set  
**Output:** VCF files with SV calls  

---

### Use Case 2: Full Benchmarking Pipeline

```yaml
# REQUIRED
fasta: '/path/to/reference.fa'
illumina_wes_bam: '/path/to/sample.bam'
benchmark_vcf: '/path/to/truth.vcf.gz'
high_confidence_targets: '/path/to/HG002_targets.bed'
gene_panel_targets: '/path/to/gene_panel.bed'
wes_utr_targets: '/path/to/wes_utr.bed'

# RECOMMENDED
skip_benchmarking: false
run_name: 'benchmarking_run'
truvari_refdist: 500
truvari_pctsize: 0.7
max_cpus: 24
max_memory: '128.GB'
```

**Purpose:** Call SVs and compare to truth set  
**Output:** VCF files + Truvari benchmarking results  

---

### Use Case 3: All Technologies (WES, WGS, PacBio, ONT)

```yaml
# REQUIRED
fasta: '/path/to/reference.fa'
illumina_wes_bam: '/path/to/wes.bam'
illumina_wgs_bam: '/path/to/wgs.bam'
pacbio_bam: '/path/to/pacbio.bam'
ont_bam: '/path/to/ont.bam'
benchmark_vcf: '/path/to/truth.vcf.gz'
high_confidence_targets: '/path/to/targets.bed'
gene_panel_targets: '/path/to/gene_panel.bed'
wes_utr_targets: '/path/to/wes_utr.bed'
tandem_repeats: '/path/to/tandem_repeats.bed'

# RECOMMENDED
skip_benchmarking: false
skip_pbsv: false
run_name: 'all_technologies'
max_cpus: 24
max_memory: '128.GB'
```

**Purpose:** Comprehensive analysis of all sequencing technologies  
**Output:** SV calls from all tools + benchmarking for each  

---

### Use Case 4: Simulation and Statistical Analysis

```yaml
# REQUIRED
fasta: '/path/to/reference.fa'
illumina_wes_bam: '/path/to/wes.bam'
benchmark_vcf: '/path/to/truth.vcf.gz'
gencode_gtf: '/path/to/gencode.v19.annotation.gtf.gz'
high_confidence_targets: '/path/to/targets.bed'
gene_panel_targets: '/path/to/gene_panel.bed'
wes_utr_targets: '/path/to/wes_utr.bed'

# RECOMMENDED
simulate_targets: true
num_simulations: 100
gather_statistics: true
skip_benchmarking: false
run_name: 'simulation_analysis'
```

**Purpose:** Simulate target regions and analyze statistical properties  
**Output:** SV calls + benchmarking + simulation results + plots  

---

### Use Case 5: Download Complete Dataset

```yaml
# REQUIRED
prepare_complete_data: true
genome: 'hs37d5'

# OPTIONAL
skip_bam_download: false
skip_reference_download: false
skip_singularity_download: false
```

**Purpose:** Download all data (BAMs, reference, benchmarks, annotations)  
**Output:** Complete dataset in `data/` directory  

---

## 🎯 Parameter Decision Tree

```
Do you have BAM files?
├─ NO → Use: prepare_complete_data: true
└─ YES → Do you have truth VCF?
    ├─ NO → Use: skip_benchmarking: true
    └─ YES → Do you have target BED files?
        ├─ NO → Create dummy BEDs OR Use: prepare_giab_resources: true
        └─ YES → Ready to run full pipeline!
```

---

## ⚙️ Tools Used per Technology

### Illumina WES
- **Manta** (germline SV caller)
- Truvari benchmarking (if enabled)

### Illumina WGS
- **Manta** (germline SV caller)
- Truvari benchmarking (if enabled)

### PacBio
- **CuteSV** (long-read SV caller)
- **PBSV** (PacBio-specific caller, unless skip_pbsv=true)
- Truvari benchmarking (if enabled)

### Oxford Nanopore (ONT)
- **CuteSV** (long-read SV caller)
- **Sniffles** (long-read SV caller)
- Truvari benchmarking (if enabled)
- Tandem repeat annotation (if tandem_repeats provided)

---

## 📝 Example Complete params.yaml

Here's a fully annotated example:

```yaml
# =====================================================
# Complete Production Parameters Example
# =====================================================

# Run identification
run_name: 'HG002_complete_analysis'
outdir: '/data/projects/sv_calling/HG002_complete'

# Reference genome (REQUIRED)
fasta: '/data/references/hs37d5.fa'

# Input BAMs (at least one REQUIRED)
illumina_wes_bam: '/data/bams/HG002_wes.bam'
illumina_wgs_bam: '/data/bams/HG002_wgs.bam'
pacbio_bam: '/data/bams/HG002_pacbio.bam'
ont_bam: '/data/bams/HG002_ont.bam'

# Benchmarking (REQUIRED for benchmarking mode)
benchmark_vcf: '/data/giab/HG002_SVs_Tier1_v0.6.vcf.gz'
skip_benchmarking: false

# Target regions (ALWAYS REQUIRED)
high_confidence_targets: '/data/giab/HG002_SVs_Tier1_v0.6.bed'
gene_panel_targets: '/data/targets/gene_panel.bed'
wes_utr_targets: '/data/targets/wes_utr_regions.bed'
tandem_repeats: '/data/annotations/tandem_repeats.bed'

# Pipeline behavior
skip_pbsv: false
simulate_targets: false
gather_statistics: true

# Truvari parameters
truvari_refdist: 500
truvari_pctsize: 0.7
truvari_pctovl: 0.0
truvari_pctseq: 0.0

# Resources
max_cpus: 32
max_memory: '256.GB'
max_time: '48.h'
```

---

## 🔍 Validation

Before running, verify your parameters:

```bash
# Check if files exist
test -f /path/to/reference.fa && echo "✓ Reference found" || echo "✗ Missing"
test -f /path/to/sample.bam && echo "✓ BAM found" || echo "✗ Missing"
test -f /path/to/sample.bam.bai && echo "✓ BAI found" || echo "✗ Missing"
test -f /path/to/targets.bed && echo "✓ Targets found" || echo "✗ Missing"
```

---

## 💡 Pro Tips

1. **Use absolute paths** - Avoid relative paths in params files
2. **Start minimal** - Test with one technology before running all
3. **Check indexes** - All BAMs need .bai, reference needs .fai, VCF needs .tbi
4. **Use dummy BEDs** - For testing without benchmarking
5. **Set resource limits** - Prevent out-of-memory errors
6. **Use -resume** - Continue from failures without rerunning completed steps

---

## 🆘 Quick Reference

```yaml
# Minimal (SV calling only)
fasta: required
*_bam: at least one required
*_targets: all three required (can be dummy)
skip_benchmarking: true

# With benchmarking
benchmark_vcf: required
*_targets: all three required (must be real)
skip_benchmarking: false

# All optional
Everything else has sensible defaults!
```
