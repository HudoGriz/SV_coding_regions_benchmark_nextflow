# Complete Data Preparation Guide

This guide explains how to use the complete data preparation workflows that replicate the functionality of the original bash scripts.

## Overview

The pipeline provides two complete data preparation workflows:
1. **GRCh37/hs37d5**: Includes Singularity containers, all BAM files, reference genome, and annotations
2. **GRCh38**: Includes BAM files, reference genome, and annotations (uses pre-downloaded containers)

These workflows automate all downloads from the original bash scripts into modular Nextflow processes.

## Quick Start

### Complete GRCh37 Setup
```bash
nextflow run main.nf -profile complete_grch37
```

### Complete GRCh38 Setup
```bash
nextflow run main.nf -profile complete_grch38
```

## What Gets Downloaded

### GRCh37/hs37d5 (`complete_grch37`)

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
- Illumina WES: `151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam`
- Illumina WGS: `HG002.hs37d5.60x.1.bam`
- PacBio HiFi: `HG002_PacBio-HiFi-Revio_20231031_48x_GRCh37.bam`
- ONT Ultralong: `HG002_GRCh37_ONT-UL_UCSC_20200508.phased.bam`

**Reference Genome:**
- hs37d5.fa (1000 Genomes Phase 2)
- hs37d5.fa.fai (index)

**GIAB Truth Sets:**
- HG002_SVs_Tier1_v0.6.vcf.gz
- HG002_SVs_Tier1_v0.6.vcf.gz.tbi
- HG002_SVs_Tier1_v0.6.bed

**Annotations:**
- human_hs37d5.trf.bed (tandem repeats)
- gencode.v19.annotation.gtf.gz

**Generated Files:**
- exome_utr_gtf.bed
- exome_utr_gtf.HG002_SVs_Tier1.bed (intersected)
- Paediatric_disorders.HG002_SVs_Tier1.bed (if panel provided)

### GRCh38 (`complete_grch38`)

**BAM Files:**
- Illumina WGS: `HG002.GRCh38.300x.bam`
- PacBio HiFi: `HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam`
- ONT Ultralong: `HG002_GRCh38_ONT-UL_UCSC_20200508.phased.bam`

**Reference Genome:**
- human_GRCh38_no_alt_analysis_set.fasta
- human_GRCh38_no_alt_analysis_set.fasta.fai

**GIAB Truth Sets:**
- GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz (native)
- GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz.tbi
- GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed
- GRCh37_HG002-T2TQ100-V1.0_stvar.vcf.gz (liftover, optional)
- GRCh37_HG002-T2TQ100-V1.0_stvar.vcf.gz.tbi
- GRCh37_HG002-T2TQ100-V1.0_stvar.benchmark.bed

**Annotations:**
- human_GRCh38_no_alt_analysis_set.trf.bed
- gencode.v49.annotation.gtf.gz

**Generated Files:**
- exome_utr_gtf_GRCh38.bed
- exome_utr_gtf.GRCh38_HG002-T2TQ100-V1.0_stvar.bed (intersected)
- Paediatric_disorders_GRCh38.GRCh38_HG002-T2TQ100-V1.0_stvar.bed (if panel provided)

## Output Directory Structure

### GRCh37
```
project_dir/
├── singularity_images/
│   ├── manta_latest.sif
│   ├── samtools_latest.sif
│   ├── cutesv_latest.sif
│   ├── pbsv_latest.sif
│   ├── sniffles_latest.sif
│   ├── bedtools_latest.sif
│   ├── truvari_modded.sif
│   └── r-env_4-4-1.sif
├── data/
│   ├── Illumina_wes/
│   │   └── bam/
│   │       ├── *.bam
│   │       └── *.bai
│   ├── Illumina_wgs/
│   │   └── bam/
│   │       ├── HG002.hs37d5.60x.1.bam
│   │       └── HG002.hs37d5.60x.1.bam.bai
│   ├── Pacbio/
│   │   └── bam/
│   │       ├── HG002_PacBio-HiFi-Revio_*.bam
│   │       └── HG002_PacBio-HiFi-Revio_*.bam.bai
│   ├── ONT/
│   │   └── bam/
│   │       ├── HG002_GRCh37_ONT-UL_*.bam
│   │       └── HG002_GRCh37_ONT-UL_*.bam.bai
│   └── references/
│       ├── human_hs37d5.fasta
│       ├── human_hs37d5.fasta.fai
│       ├── HG002_SVs_Tier1_v0.6.vcf.gz
│       ├── HG002_SVs_Tier1_v0.6.vcf.gz.tbi
│       ├── HG002_SVs_Tier1_v0.6.bed
│       ├── human_hs37d5.trf.bed
│       ├── gencode.v19.annotation.gtf.gz
│       ├── exome_utr_gtf.bed
│       └── exome_utr_gtf.HG002_SVs_Tier1.bed
```

### GRCh38
```
project_dir/
├── data/
│   ├── Illumina_wgs/
│   │   └── bam_GRCh38/
│   │       ├── HG002.GRCh38.300x.bam
│   │       └── HG002.GRCh38.300x.bam.bai
│   ├── Pacbio/
│   │   └── bam_GRCh38/
│   │       ├── HG002_PacBio-HiFi-Revio_*_GRCh38-GIABv3.bam
│   │       └── HG002_PacBio-HiFi-Revio_*_GRCh38-GIABv3.bam.bai
│   ├── ONT/
│   │   └── bam_GRCh38/
│   │       ├── HG002_GRCh38_ONT-UL_*.bam
│   │       └── HG002_GRCh38_ONT-UL_*.bam.bai
│   └── references/
│       ├── human_GRCh38_no_alt_analysis_set.fasta
│       ├── human_GRCh38_no_alt_analysis_set.fasta.fai
│       ├── GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz
│       ├── GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz.tbi
│       ├── GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed
│       ├── GRCh37_HG002-T2TQ100-V1.0_stvar.vcf.gz
│       ├── GRCh37_HG002-T2TQ100-V1.0_stvar.vcf.gz.tbi
│       ├── GRCh37_HG002-T2TQ100-V1.0_stvar.benchmark.bed
│       ├── human_GRCh38_no_alt_analysis_set.trf.bed
│       ├── gencode.v49.annotation.gtf.gz
│       ├── exome_utr_gtf_GRCh38.bed
│       └── exome_utr_gtf.GRCh38_HG002-T2TQ100-V1.0_stvar.bed
```

## Advanced Usage

### Skip Already Downloaded Files

If you already have some files, you can skip their download:

```bash
# Skip Singularity containers (GRCh37 only)
nextflow run main.nf -profile complete_grch37 --skip_singularity_download

# Skip BAM files
nextflow run main.nf -profile complete_grch37 --skip_bam_download

# Skip reference genome
nextflow run main.nf -profile complete_grch37 --skip_reference_download

# Skip multiple
nextflow run main.nf -profile complete_grch37 \
    --skip_singularity_download \
    --skip_bam_download
```

### Skip GRCh37 Liftover (GRCh38 only)

```bash
nextflow run main.nf -profile complete_grch38 --download_grch37_liftover false
```

### Add Gene Panel Intersections

If you have custom gene panel BED files:

```bash
# GRCh37
nextflow run main.nf -profile complete_grch37 \
    --paediatric_disorders_bed /path/to/Paediatric_disorders.bed

# GRCh38
nextflow run main.nf -profile complete_grch38 \
    --paediatric_disorders_bed_grch38 /path/to/Paediatric_disorders_GRCh38.bed
```

### Custom Project Directory

```bash
nextflow run main.nf -profile complete_grch37 \
    --project_dir /path/to/custom/location
```

## Parameters Reference

| Parameter | Default | Description |
|-----------|---------|-------------|
| `prepare_complete_data` | `false` | Enable complete data preparation |
| `genome` | `'hs37d5'` | Genome build (`'hs37d5'` or `'GRCh38'`) |
| `project_dir` | `${projectDir}` | Base directory for outputs |
| `skip_singularity_download` | `false` | Skip Singularity container downloads |
| `skip_bam_download` | `false` | Skip BAM file downloads |
| `skip_reference_download` | `false` | Skip reference genome download |
| `download_grch37_liftover` | `true` | Download GRCh37 liftover (GRCh38 only) |
| `paediatric_disorders_bed` | `null` | Custom gene panel BED (GRCh37) |
| `paediatric_disorders_bed_grch38` | `null` | Custom gene panel BED (GRCh38) |
| `r_container` | `${projectDir}/singularity_images/r-env_4-4-1.sif` | R container path |

## Resource Requirements

### GRCh37 Complete
- **Download Size**: ~200 GB (BAMs are large)
- **Storage Required**: ~250 GB
- **Runtime**: 4-8 hours (depends on network speed)
- **Peak Memory**: 8 GB

### GRCh38 Complete
- **Download Size**: ~180 GB
- **Storage Required**: ~230 GB
- **Runtime**: 4-8 hours
- **Peak Memory**: 8 GB

## Data Sources

### Singularity Containers
- Docker Hub (manta, samtools, cutesv, sniffles, bedtools)
- Quay.io (PacBio pbsv)
- Sylabs Cloud (truvari_modded, r-env)

### BAM Files
- NCBI FTP GIAB Repository
- Various sequencing technologies and depths

### Reference Genomes
- 1000 Genomes (hs37d5)
- NCBI (GRCh38 no_alt analysis set)

### Truth Sets
- GIAB v0.6 (GRCh37)
- GIAB T2TQ100-V1.0 (GRCh38 and liftover)

### Annotations
- PacBio pbsv repository (tandem repeats)
- GENCODE (gene annotations)

## Comparison with Bash Scripts

| Bash Script Feature | Nextflow Implementation | Benefits |
|---------------------|------------------------|----------|
| Sequential downloads | Parallel processes | Faster execution |
| Manual error checking | Automatic retries | More robust |
| Hard-coded paths | Parameterized | Flexible |
| Single script | Modular processes | Reusable components |
| No resume capability | `-resume` support | Efficient restarts |
| Manual skip logic | Skip parameters | Easy customization |

## Troubleshooting

### Download Failures

Nextflow automatically retries failed downloads. To resume after errors:
```bash
nextflow run main.nf -profile complete_grch37 -resume
```

### Disk Space Errors

Check available space before starting:
```bash
df -h $PWD
```

Ensure at least 250 GB free for complete downloads.

### Network Timeouts

For slow networks, increase retry attempts in `nextflow.config`:
```groovy
withName: 'DOWNLOAD_.*' {
    errorStrategy = 'retry'
    maxRetries = 5  // Increase from 3
}
```

### Container Pull Failures

If Singularity pulls fail, manually pull problematic containers:
```bash
cd singularity_images
singularity pull manta_latest.sif docker://dceoy/manta:latest
```

Then skip container downloads:
```bash
nextflow run main.nf -profile complete_grch37 --skip_singularity_download
```

### Missing R Packages

If CREATE_EXOME_UTR_BED fails, verify R container has required packages:
```bash
singularity exec singularity_images/r-env_4-4-1.sif R -e 'installed.packages()[,1]'
```

## Integration Example

After complete data preparation, run the SV calling pipeline:

```bash
# Step 1: Prepare all data
nextflow run main.nf -profile complete_grch37

# Step 2: Run SV calling (example with Illumina WES)
nextflow run main.nf \
    --illumina_wes_bam data/Illumina_wes/bam/*.bam \
    --fasta data/references/human_hs37d5.fasta \
    --benchmark_vcf data/references/HG002_SVs_Tier1_v0.6.vcf.gz \
    --high_confidence_targets data/references/HG002_SVs_Tier1_v0.6.bed \
    --wes_utr_targets data/references/exome_utr_gtf.bed \
    --outdir results/sv_calling
```

## Notes

- The workflows are idempotent - safe to run multiple times
- Use `-resume` to continue interrupted downloads
- Skip parameters allow incremental setup
- GRCh37 includes WES BAMs, GRCh38 does not (comment out in script)
- All intermediate files are preserved for reproducibility
- Downloaded files persist outside `work/` directory

## Next Steps

After data preparation:
1. Review downloaded files in `data/` directory
2. Verify file integrity (check file sizes, MD5 sums if available)
3. Run SV calling pipeline with prepared data
4. Use prepared truth sets for benchmarking with Truvari

## Support

For issues or questions:
- Check workflow logs in `.nextflow.log`
- Review process logs in `work/` directories
- Verify URLs are still valid (NCBI/GENCODE may update)
- Ensure sufficient disk space and network bandwidth
