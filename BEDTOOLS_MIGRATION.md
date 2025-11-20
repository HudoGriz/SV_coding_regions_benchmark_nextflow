# Bedtools Module Migration to nf-core

## Summary
Successfully replaced local `bedtools_intersect` module with nf-core `bedtools/intersect` module across all workflow files.

## Changes Made

### 1. Updated Workflow Files

#### workflows/preparation/prepare_data_complete_grch37.nf
- ✅ Updated includes from `../../modules/local/bedtools_intersect` to `../../modules/nf-core/bedtools/intersect/main`
- ✅ Updated process calls to match nf-core signature:
  - Added proper channel construction with `.map` operations
  - Added empty chrom_sizes channel: `[[id: 'null'], []]`
- ✅ Updated output references from `.out.bed` to `.out.intersect`

#### workflows/preparation/prepare_data_complete_grch38.nf
- ✅ Updated includes from `../../modules/local/bedtools_intersect` to `../../modules/nf-core/bedtools/intersect/main`
- ✅ Updated process calls to match nf-core signature
- ✅ Updated output references from `.out.bed` to `.out.intersect`

#### workflows/prepare_data_grch37.nf
- ✅ Updated process calls to properly construct channels
- ✅ Updated output references from `.out.bed` to `.out.intersect`

#### workflows/prepare_data_grch38.nf
- ✅ Updated process calls to properly construct channels
- ✅ Updated output references from `.out.bed` to `.out.intersect`

### 2. Configuration Updates

#### nextflow.config
- ✅ Added `ext.suffix = 'bed'` to ensure .bed file extension
- ✅ Added publishDir configuration for BEDTOOLS_INTERSECT processes
- ✅ Updated pattern matching to include all INTERSECT_* process names

## Key Differences Between Local and nf-core Modules

### Local Module (bedtools_intersect.nf)
```groovy
input:
tuple val(meta), path(bed_a), path(bed_b)

output:
tuple val(meta), path("*.bed"), emit: bed

container 'file:///singularity_images/bedtools_latest.sif'
```

### nf-core Module (bedtools/intersect)
```groovy
input:
tuple val(meta), path(intervals1), path(intervals2)
tuple val(meta2), path(chrom_sizes)

output:
tuple val(meta), path("*.${extension}"), emit: intersect
path  "versions.yml"                   , emit: versions

container (auto-selected from Galaxy depot or biocontainers)
```

## Benefits of Migration

1. **No Local Singularity Containers**: nf-core modules automatically pull containers from public registries
2. **Version Tracking**: Built-in version tracking in `versions.yml`
3. **Standardization**: Follows nf-core best practices
4. **Flexibility**: `ext.suffix` allows easy output format control
5. **Maintainability**: Updates from nf-core community

## Testing Recommendations

Before running the full pipeline, test with:
```bash
# Test data preparation workflow
nextflow run main.nf -profile complete_grch37 --prepare_complete_data true --skip_bam_download true
```

## Next Steps

The local module file can now be safely removed:
- `modules/local/bedtools_intersect.nf` (no longer referenced)

## Migration Status
- ✅ bedtools - **COMPLETED**
- ⏳ gunzip - In progress
- ⏳ tabix - Pending
- ⏳ samtools - Pending
