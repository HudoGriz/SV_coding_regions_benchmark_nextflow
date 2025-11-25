# Requirements Checklist for Running the Pipeline

Before running the pipeline, make sure you have all required files and their indexes.

## ✅ Required Files for Minimal Run (SV Calling Only)

### 1. Reference Genome
- [ ] **FASTA file**: `data/reference/hs37d5.fa` (or GRCh37/GRCh38)
- [ ] **FASTA index**: `data/reference/hs37d5.fa.fai` (auto-created if missing)

**How to create index:**
```bash
samtools faidx data/reference/hs37d5.fa
```

---

### 2. Input BAM Files (at least ONE required)

#### Illumina WES
- [ ] **BAM file**: `data/Illumina_wes/bam/*.bam`
- [ ] **BAM index**: `data/Illumina_wes/bam/*.bam.bai`

#### Illumina WGS
- [ ] **BAM file**: `data/Illumina_wgs/bam/*.bam`
- [ ] **BAM index**: `data/Illumina_wgs/bam/*.bam.bai`

#### PacBio
- [ ] **BAM file**: `data/Pacbio/bam/*.bam`
- [ ] **BAM index**: `data/Pacbio/bam/*.bam.bai`

#### Oxford Nanopore (ONT)
- [ ] **BAM file**: `data/ONT/bam/*.bam`
- [ ] **BAM index**: `data/ONT/bam/*.bam.bai`

**How to create BAM indexes:**
```bash
samtools index data/Illumina_wes/bam/sample.bam
samtools index data/Illumina_wgs/bam/sample.bam
samtools index data/Pacbio/bam/sample.bam
samtools index data/ONT/bam/sample.bam
```

---

### 3. Target BED Files (ALWAYS REQUIRED - even without benchmarking!)

**⚠️ CRITICAL:** These are required by the workflow even if `skip_benchmarking: true`

- [ ] **High confidence targets**: `data/targets/HG002_SVs_Tier1_v0.6.bed`
- [ ] **Gene panel targets**: `data/targets/gene_panel.bed`
- [ ] **WES + UTR targets**: `data/targets/wes_utr_regions.bed`

**Option A: Use actual target files**
```bash
# If you have GIAB/benchmark data
ls -lh data/targets/*.bed
```

**Option B: Create dummy BED files for testing**
```bash
# Use the helper script
chmod +x create_dummy_targets.sh
./create_dummy_targets.sh

# Or manually create dummy files
mkdir -p data/targets
echo -e "chr1\t1\t249250621" > data/targets/dummy_high_confidence.bed
echo -e "chr1\t1\t249250621" > data/targets/dummy_gene_panel.bed
echo -e "chr1\t1\t249250621" > data/targets/dummy_wes_utr.bed
```

---

## ✅ Additional Files for Benchmarking

If you want to enable benchmarking (`skip_benchmarking: false`):

### 4. Truth/Benchmark VCF
- [ ] **Benchmark VCF**: `data/benchmark/HG002_SVs_Tier1_v0.6.vcf.gz`
- [ ] **VCF index**: `data/benchmark/HG002_SVs_Tier1_v0.6.vcf.gz.tbi`

**How to create VCF index:**
```bash
tabix -p vcf data/benchmark/HG002_SVs_Tier1_v0.6.vcf.gz
```

---

## ✅ Optional Files

### 5. Tandem Repeats (for ONT/Sniffles)
- [ ] **Tandem repeats BED**: `data/targets/tandem_repeats.bed`

Only needed if running ONT data with Sniffles and want to annotate tandem repeats.

### 6. GENCODE GTF (for simulation mode)
- [ ] **GENCODE annotation**: `data/annotations/gencode.v19.annotation.gtf.gz`

Only needed if `simulate_targets: true`.

---

## 📋 Quick Validation Script

Run this to check if all your files exist:

```bash
#!/bin/bash

PROJECT_DIR="/home/45483vrhovsek/genome_in_bottle/SV_coding_regions_benchmark_nextflow"

echo "Checking required files..."
echo ""

# Reference
echo "1. Reference genome:"
[ -f "$PROJECT_DIR/data/reference/hs37d5.fa" ] && echo "  ✓ FASTA found" || echo "  ✗ FASTA NOT FOUND"
[ -f "$PROJECT_DIR/data/reference/hs37d5.fa.fai" ] && echo "  ✓ FAI index found" || echo "  ⚠ FAI index missing (will be created)"
echo ""

# BAM files
echo "2. BAM files:"
if [ -f "$PROJECT_DIR/data/Illumina_wes/bam/"*.bam ]; then
    echo "  ✓ Illumina WES BAM found"
    [ -f "$PROJECT_DIR/data/Illumina_wes/bam/"*.bam.bai ] && echo "  ✓ BAI index found" || echo "  ✗ BAI index NOT FOUND"
else
    echo "  ✗ Illumina WES BAM not found"
fi
echo ""

# Target BED files
echo "3. Target BED files (REQUIRED):"
[ -f "$PROJECT_DIR/data/targets/HG002_SVs_Tier1_v0.6.bed" ] && echo "  ✓ High confidence targets" || echo "  ✗ High confidence targets NOT FOUND"
[ -f "$PROJECT_DIR/data/targets/gene_panel.bed" ] && echo "  ✓ Gene panel targets" || echo "  ✗ Gene panel targets NOT FOUND"
[ -f "$PROJECT_DIR/data/targets/wes_utr_regions.bed" ] && echo "  ✓ WES+UTR targets" || echo "  ✗ WES+UTR targets NOT FOUND"
echo ""

# Benchmark VCF (optional)
echo "4. Benchmark VCF (optional, for benchmarking):"
if [ -f "$PROJECT_DIR/data/benchmark/HG002_SVs_Tier1_v0.6.vcf.gz" ]; then
    echo "  ✓ Benchmark VCF found"
    [ -f "$PROJECT_DIR/data/benchmark/HG002_SVs_Tier1_v0.6.vcf.gz.tbi" ] && echo "  ✓ TBI index found" || echo "  ✗ TBI index NOT FOUND"
else
    echo "  ⚠ Benchmark VCF not found (OK if skip_benchmarking=true)"
fi
echo ""
```

---

## 🚦 Minimum Requirements Summary

### To run SV calling only (no benchmarking):

```yaml
# REQUIRED FILES:
fasta: reference.fa + .fai
illumina_wes_bam: sample.bam + .bai  # at least one BAM
high_confidence_targets: dummy.bed   # can be dummy
gene_panel_targets: dummy.bed        # can be dummy
wes_utr_targets: dummy.bed          # can be dummy

# SETTINGS:
skip_benchmarking: true
```

### To run with benchmarking:

```yaml
# ALL ABOVE PLUS:
benchmark_vcf: truth.vcf.gz + .tbi
high_confidence_targets: real.bed    # must be real
gene_panel_targets: real.bed         # must be real
wes_utr_targets: real.bed           # must be real

# SETTINGS:
skip_benchmarking: false
```

---

## 🔍 Common Issues

### Issue: "File not found" error
**Solution:** Check that all paths in params file are absolute and correct

### Issue: "Index not found" error
**Solution:** Create missing indexes:
```bash
samtools faidx reference.fa
samtools index sample.bam
tabix -p vcf truth.vcf.gz
```

### Issue: "Cannot access null object"
**Solution:** Make sure all three target BED files are specified (can be dummy files)

### Issue: Pipeline fails at PREPARE_REFERENCES
**Solution:** Check that your reference FASTA exists and is readable

---

## 📁 Expected Directory Structure

```
/home/45483vrhovsek/genome_in_bottle/SV_coding_regions_benchmark_nextflow/
├── main.nf
├── params_minimal_example.yaml
├── params_production.yaml
├── data/
│   ├── reference/
│   │   ├── hs37d5.fa
│   │   └── hs37d5.fa.fai
│   ├── Illumina_wes/bam/
│   │   ├── sample.bam
│   │   └── sample.bam.bai
│   ├── Illumina_wgs/bam/
│   │   ├── sample.bam
│   │   └── sample.bam.bai
│   ├── Pacbio/bam/
│   │   ├── sample.bam
│   │   └── sample.bam.bai
│   ├── ONT/bam/
│   │   ├── sample.bam
│   │   └── sample.bam.bai
│   ├── targets/
│   │   ├── HG002_SVs_Tier1_v0.6.bed
│   │   ├── gene_panel.bed
│   │   └── wes_utr_regions.bed
│   └── benchmark/
│       ├── HG002_SVs_Tier1_v0.6.vcf.gz
│       └── HG002_SVs_Tier1_v0.6.vcf.gz.tbi
└── work/  (created by Nextflow)
```

---

## ✅ Ready to Run Checklist

Before executing `nextflow run main.nf -params-file params_minimal_example.yaml`:

- [ ] Reference FASTA exists (index will be auto-created)
- [ ] At least one BAM file exists with .bai index
- [ ] All three target BED files exist (can be dummy files)
- [ ] All paths in params file are correct and absolute
- [ ] Singularity or Docker is installed and working
- [ ] Sufficient disk space for outputs (~50GB per technology)
- [ ] Sufficient memory (8GB minimum, 32GB recommended)

Once all checked, you're ready to run! 🚀
