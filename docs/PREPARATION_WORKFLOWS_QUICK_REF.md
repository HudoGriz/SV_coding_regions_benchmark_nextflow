# Data Preparation Workflows - Quick Reference

## Available Workflows

The pipeline provides three levels of data preparation:

| Workflow | What It Does | Use When |
|----------|--------------|----------|
| **GIAB Resources Only** | Downloads truth sets & annotations only | You have BAMs/reference already |
| **Complete GRCh37** | Everything including Singularity containers | Fresh setup for GRCh37 |
| **Complete GRCh38** | Everything (uses pre-downloaded containers) | Fresh setup for GRCh38 |

## One-Command Usage

### Minimal: GIAB Resources Only

```bash
# GRCh37 truth sets & annotations
nextflow run main.nf -profile giab_grch37

# GRCh38 truth sets & annotations
nextflow run main.nf -profile giab_grch38
```

**Downloads**: Truth VCF, high-confidence BED, tandem repeats, GENCODE, exome+UTR BED  
**Time**: 10-30 minutes  
**Size**: ~50-75 MB

### Complete: Everything You Need

```bash
# Complete GRCh37 setup (includes Singularity containers)
nextflow run main.nf -profile complete_grch37

# Complete GRCh38 setup (uses pre-downloaded containers)
nextflow run main.nf -profile complete_grch38
```

**GRCh37 Downloads**: Containers + BAMs + Reference + Truth sets + Annotations  
**GRCh38 Downloads**: BAMs + Reference + Truth sets + Annotations  
**Time**: 3-8 hours  
**Size**: 200-250 GB

## What Each Profile Downloads

### `giab_grch37` (Minimal)
- ✅ HG002_SVs_Tier1_v0.6.vcf.gz + index
- ✅ HG002_SVs_Tier1_v0.6.bed
- ✅ human_hs37d5.trf.bed
- ✅ gencode.v19.annotation.gtf.gz
- ✅ exome_utr_gtf.bed (generated)

### `giab_grch38` (Minimal)
- ✅ GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz + index
- ✅ GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed
- ✅ human_GRCh38_no_alt_analysis_set.trf.bed
- ✅ gencode.v49.annotation.gtf.gz
- ✅ exome_utr_gtf_GRCh38.bed (generated)

### `complete_grch37` (Full Setup)
**Singularity Containers:**
- manta_latest.sif
- samtools_latest.sif
- cutesv_latest.sif
- pbsv_latest.sif
- sniffles_latest.sif
- bedtools_latest.sif
- truvari_modded.sif
- r-env_4-4-1.sif

**BAM Files:**
- Illumina WES: ~15 GB
- Illumina WGS: ~120 GB
- PacBio HiFi: ~50 GB
- ONT Ultralong: ~80 GB

**Reference + Everything from giab_grch37**

### `complete_grch38` (Full Setup)
**BAM Files:**
- Illumina WGS: ~180 GB
- PacBio HiFi: ~50 GB
- ONT Ultralong: ~70 GB

**Reference:**
- human_GRCh38_no_alt_analysis_set.fasta

**Everything from giab_grch38** PLUS:
- GRCh37_HG002-T2TQ100-V1.0_stvar.vcf.gz (liftover)
- GRCh37_HG002-T2TQ100-V1.0_stvar.benchmark.bed

## Advanced Options

### Skip Already Downloaded Files

```bash
# Skip Singularity containers (GRCh37 only)
nextflow run main.nf -profile complete_grch37 --skip_singularity_download

# Skip BAM files
nextflow run main.nf -profile complete_grch37 --skip_bam_download

# Skip reference genome
nextflow run main.nf -profile complete_grch37 --skip_reference_download

# Combine multiple
nextflow run main.nf -profile complete_grch37 \
    --skip_singularity_download \
    --skip_bam_download
```

### Resume Interrupted Downloads

```bash
# If download was interrupted, resume:
nextflow run main.nf -profile complete_grch37 -resume
```

### Custom Output Location

```bash
nextflow run main.nf -profile complete_grch37 \
    --project_dir /path/to/custom/location
```

### Add Gene Panel Intersections

```bash
# GRCh37
nextflow run main.nf -profile complete_grch37 \
    --paediatric_disorders_bed /path/to/panel.bed

# GRCh38
nextflow run main.nf -profile complete_grch38 \
    --paediatric_disorders_bed_grch38 /path/to/panel_grch38.bed
```

## Output Locations

### After `giab_grch37` or `complete_grch37`
```
project_dir/data/references/
├── HG002_SVs_Tier1_v0.6.vcf.gz
├── HG002_SVs_Tier1_v0.6.bed
├── human_hs37d5.trf.bed
├── gencode.v19.annotation.gtf.gz
└── exome_utr_gtf.bed
```

### After `giab_grch38` or `complete_grch38`
```
project_dir/data/references/
├── GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz
├── GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed
├── GRCh37_HG002-T2TQ100-V1.0_stvar.vcf.gz (liftover)
├── GRCh37_HG002-T2TQ100-V1.0_stvar.benchmark.bed
├── human_GRCh38_no_alt_analysis_set.trf.bed
├── gencode.v49.annotation.gtf.gz
└── exome_utr_gtf_GRCh38.bed
```

### Additional Outputs from `complete_*` Profiles
```
project_dir/
├── singularity_images/           (GRCh37 only)
│   └── *.sif
└── data/
    ├── Illumina_wes/bam/          (GRCh37 only)
    ├── Illumina_wgs/bam/          (GRCh37)
    ├── Illumina_wgs/bam_GRCh38/   (GRCh38)
    ├── Pacbio/bam/                (GRCh37)
    ├── Pacbio/bam_GRCh38/         (GRCh38)
    ├── ONT/bam/                   (GRCh37)
    ├── ONT/bam_GRCh38/            (GRCh38)
    └── references/
        ├── human_hs37d5.fasta     (GRCh37)
        └── human_GRCh38_no_alt_analysis_set.fasta (GRCh38)
```

## Workflow Selection Guide

### Choose `giab_grch37` when:
- ✅ You already have BAM files
- ✅ You already have reference genome
- ✅ You only need truth sets for benchmarking
- ✅ Quick setup (10-30 min)

### Choose `complete_grch37` when:
- ✅ Fresh installation
- ✅ Need Singularity containers
- ✅ Need all HG002 BAM files
- ✅ Need reference genome
- ✅ Full automated setup

### Choose `giab_grch38` when:
- ✅ Have GRCh38 BAMs and reference
- ✅ Only need GRCh38 truth sets
- ✅ Quick setup (10-30 min)

### Choose `complete_grch38` when:
- ✅ Fresh GRCh38 installation
- ✅ Need all HG002 GRCh38 BAMs
- ✅ Need GRCh38 reference
- ✅ Want both native and liftover truth sets

## Resource Requirements

| Workflow | Download Size | Final Size | Time | Memory | CPUs |
|----------|--------------|------------|------|--------|------|
| `giab_grch37` | ~50 MB | ~500 MB | 10-30 min | 8 GB | 1-2 |
| `giab_grch38` | ~75 MB | ~800 MB | 10-30 min | 8 GB | 1-2 |
| `complete_grch37` | ~200 GB | ~250 GB | 3-8 hrs | 8 GB | 4-8 |
| `complete_grch38` | ~180 GB | ~230 GB | 3-8 hrs | 8 GB | 4-8 |

## Common Tasks

### Check What Would Be Downloaded
```bash
# Dry run (doesn't actually download)
nextflow run main.nf -profile complete_grch37 -stub-run
```

### Resume After Network Failure
```bash
nextflow run main.nf -profile complete_grch37 -resume
```

### View Execution Report
```bash
nextflow run main.nf -profile complete_grch37 -with-report report.html
open report.html
```

### Clean Up Work Directory
```bash
# After successful completion
nextflow clean -f
```

## Troubleshooting Quick Fixes

### Download Failed
```bash
# Resume with automatic retry
nextflow run main.nf -profile complete_grch37 -resume
```

### Out of Disk Space
```bash
# Check available space
df -h .

# Move work directory to larger disk
nextflow run main.nf -profile complete_grch37 \
    -work-dir /path/to/large/disk/work
```

### Container Pull Failed
```bash
# Skip container download, use existing
nextflow run main.nf -profile complete_grch37 \
    --skip_singularity_download
```

### Specific URL Failed
```bash
# Check .nextflow.log for details
grep ERROR .nextflow.log

# Resume after fixing (automatic retry)
nextflow run main.nf -profile complete_grch37 -resume
```

## Next Steps After Preparation

### Using Prepared GRCh37 Data
```bash
nextflow run main.nf \
    --illumina_wes_bam data/Illumina_wes/bam/*.bam \
    --fasta data/references/human_hs37d5.fasta \
    --benchmark_vcf data/references/HG002_SVs_Tier1_v0.6.vcf.gz \
    --high_confidence_targets data/references/HG002_SVs_Tier1_v0.6.bed \
    --wes_utr_targets data/references/exome_utr_gtf.bed
```

### Using Prepared GRCh38 Data
```bash
nextflow run main.nf \
    --illumina_wgs_bam data/Illumina_wgs/bam_GRCh38/*.bam \
    --fasta data/references/human_GRCh38_no_alt_analysis_set.fasta \
    --benchmark_vcf data/references/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz \
    --high_confidence_targets data/references/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed \
    --wes_utr_targets data/references/exome_utr_gtf_GRCh38.bed
```

## Documentation Links

- **User Guide**: `docs/guides/complete_data_preparation.md`
- **GIAB Resources Guide**: `docs/guides/prepare_giab_resources.md`
- **Migration Guide**: `docs/BASH_TO_NEXTFLOW_MIGRATION.md`
- **Implementation Summary**: `COMPLETE_PREPARATION_IMPLEMENTATION.md`

## Summary Comparison

| Feature | giab_* | complete_* |
|---------|--------|------------|
| **Truth Sets** | ✅ | ✅ |
| **Annotations** | ✅ | ✅ |
| **Target BEDs** | ✅ | ✅ |
| **BAM Files** | ❌ | ✅ |
| **Reference** | ❌ | ✅ |
| **Containers** | ❌ | ✅ (GRCh37 only) |
| **Time** | 10-30 min | 3-8 hrs |
| **Size** | 50-75 MB | 200-250 GB |
| **Use Case** | Have data, need benchmarks | Full automated setup |
