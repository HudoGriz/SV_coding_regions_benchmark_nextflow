# PREPARE_GIAB_RESOURCES Workflow

## Overview
This workflow downloads and prepares Genome in a Bottle (GIAB) HG002 resources for variant calling validation. It adapts to the reference genome version (GRCh37 or GRCh38) and downloads the appropriate truth sets, annotations, and creates target region BED files.

## Workflow Logic

```
START → Check genome version → Download truth sets → Download annotations → Create target BED → END
```

### Genome-Specific Paths

**GRCh37/hs37d5:**
- Truth Set: v0.6 (legacy format)
- Tandem Repeats: hs37d5
- GENCODE: v19

**GRCh38:**
- Truth Set: T2TQ100-V1.0
- Tandem Repeats: GRCh38
- GENCODE: v49

## Inputs

| Input | Type | Description |
|-------|------|-------------|
| `params.genome` | String | Reference genome identifier (`hs37d5` or `GRCh38`) |
| `params.project_dir` | Path | Project root directory |
| `params.r_container` | Path | Path to Singularity R container |

## Outputs

| Output | Type | Description |
|--------|------|-------------|
| `giab_vcf` | Tuple | Truth set VCF + index |
| `giab_bed` | Path | High-confidence regions BED |
| `tandem_repeats` | Path | Tandem repeat regions |
| `gencode_gtf` | Path | Gene annotations |
| `exome_utr_bed` | Path | Exome + UTR target regions |

## Output Directory Structure

```
${params.project_dir}/data/HG002_references/
├── HG002_SVs_Tier1_v0.6.vcf.gz          # GRCh37 truth VCF
├── HG002_SVs_Tier1_v0.6.vcf.gz.tbi
├── HG002_SVs_Tier1_v0.6.bed
└── human_hs37d5.trf.bed

${params.project_dir}/data/references/
├── gencode.v19.annotation.gtf.gz         # GRCh37
└── exome_utr_gtf.bed
```

Or for GRCh38:

```
${params.project_dir}/data/HG002_references/
├── GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz
├── GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz.tbi
├── GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed
└── human_GRCh38_no_alt_analysis_set.trf.bed

${params.project_dir}/data/references/
├── gencode.v49.annotation.gtf.gz
└── exome_utr_gtf_GRCh38.bed
```

## Processes

### DOWNLOAD_GIAB_TRUTH_SET
Downloads GIAB HG002 v0.6 truth set (GRCh37/hs37d5).

### DOWNLOAD_GIAB_TRUTH_SET_GRCH38
Downloads GIAB HG002 T2TQ100-V1.0 truth set (GRCh38 native).

### DOWNLOAD_GIAB_TRUTH_SET_GRCH37_LIFTOVER
Downloads GIAB HG002 T2TQ100-V1.0 truth set lifted over to GRCh37.

### DOWNLOAD_TANDEM_REPEATS / DOWNLOAD_TANDEM_REPEATS_GRCH38
Downloads tandem repeat annotations for the appropriate genome build.

### DOWNLOAD_GENCODE_GTF / DOWNLOAD_GENCODE_GTF_GRCH38
Downloads GENCODE gene annotations for the appropriate genome build.

### CREATE_EXOME_UTR_BED
Creates exome + UTR target regions from GENCODE GTF using R script.

## Usage

```nextflow
include { PREPARE_GIAB_RESOURCES } from './workflows/prepare_giab_resources'

workflow {
    PREPARE_GIAB_RESOURCES()
}
```

## Dependencies

- `wget` for file downloads
- `tabix` for VCF indexing (GRCh37 only; GRCh38 downloads pre-indexed)
- R environment with required packages (for BED creation)
- Singularity container with R

## Parameters

```groovy
params.genome = 'hs37d5'  // or 'GRCh38'
params.project_dir = '/path/to/project'
params.r_container = '/path/to/r_container.sif'
```

## Notes

- The workflow automatically detects GRCh37 vs GRCh38 and routes to appropriate processes
- All downloads include retry logic for network failures
- Output directories are created if they don't exist
- For GRCh37, you can use either the v0.6 truth set or the newer T2TQ100-V1.0 lifted-over version

## Source Files

- Main workflow: `workflows/prepare_giab_resources.nf`
- Download truth sets: `modules/local/download_truth_set.nf`
- Download annotations: `modules/local/download_annotations.nf`
- Create target BEDs: `modules/local/create_target_beds.nf`
