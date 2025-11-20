# Quick Reference: GIAB Resources Preparation

## One-Command Usage

### GRCh37/hs37d5
```bash
nextflow run main.nf -profile giab_grch37
```

### GRCh38
```bash
nextflow run main.nf -profile giab_grch38
```

## What Gets Downloaded

| Resource | GRCh37 | GRCh38 |
|----------|--------|--------|
| **Truth VCF** | v0.6 | T2TQ100-V1.0 |
| **High-conf BED** | Tier1 v0.6 | T2TQ100 benchmark |
| **Tandem Repeats** | hs37d5 | GRCh38 |
| **GENCODE** | v19 | v49 |
| **Exome+UTR BED** | Generated from v19 | Generated from v49 |

## Output Locations

```
project_dir/
├── data/
│   ├── HG002_references/
│   │   ├── *truth_set.vcf.gz
│   │   ├── *truth_set.vcf.gz.tbi
│   │   ├── *benchmark.bed
│   │   └── *tandem_repeats.bed
│   └── references/
│       ├── gencode.*.gtf.gz
│       └── exome_utr_*.bed
```

## Using Prepared Resources

After preparation, reference the files in your pipeline:

```bash
nextflow run main.nf \
    --illumina_wes_bam sample.bam \
    --benchmark_vcf data/HG002_references/[truth_vcf] \
    --high_confidence_targets data/HG002_references/[benchmark_bed] \
    --wes_utr_targets data/references/[exome_utr_bed] \
    --tandem_repeats data/HG002_references/[tandem_repeats]
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--prepare_giab_resources` | `false` | Enable preparation |
| `--genome` | `'hs37d5'` | `'hs37d5'` or `'GRCh38'` |
| `--project_dir` | `${projectDir}` | Output base directory |
| `--r_container` | `${projectDir}/singularity_images/r_container.sif` | R container path |

## Requirements

- Internet connection (NCBI FTP, GENCODE, GitHub)
- Singularity with R container (rtracklayer, GenomicFeatures)
- ~500MB free space (GRCh37) or ~800MB (GRCh38)

## Troubleshooting

### Downloads fail
```bash
# Retry with resume
nextflow run main.nf -profile giab_grch37 -resume
```

### R script fails
```bash
# Check container
singularity exec r_container.sif R -e 'installed.packages()' | grep -E "rtracklayer|GenomicFeatures"
```

### Permission errors
```bash
# Ensure directory is writable
chmod -R u+w data/
```

## Workflow Steps

1. ✅ Detect genome version
2. ✅ Download GIAB truth set
3. ✅ Download tandem repeat annotations
4. ✅ Download GENCODE GTF
5. ✅ Generate exome+UTR BED from GTF
6. ✅ Publish all files to output directories

## Files Created

### GRCh37
- `HG002_SVs_Tier1_v0.6.vcf.gz` + `.tbi`
- `HG002_SVs_Tier1_v0.6.bed`
- `human_hs37d5.trf.bed`
- `gencode.v19.annotation.gtf.gz`
- `exome_utr_gtf.bed`

### GRCh38
- `GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz` + `.tbi`
- `GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed`
- `human_GRCh38_no_alt_analysis_set.trf.bed`
- `gencode.v49.annotation.gtf.gz`
- `exome_utr_gtf_GRCh38.bed`

## Data Sources

- **GIAB**: ftp-trace.ncbi.nlm.nih.gov
- **GENCODE**: ftp.ebi.ac.uk/pub/databases/gencode
- **Tandem Repeats**: github.com/hall-lab/speedseq

## Full Documentation

- User Guide: `docs/guides/prepare_giab_resources.md`
- Technical Docs: `docs/workflows/prepare_giab_resources.md`
- Implementation: `GIAB_RESOURCES_IMPLEMENTATION.md`
