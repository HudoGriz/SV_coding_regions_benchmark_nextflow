# Migration: SAMTOOLS_FAIDX and TABIX_TABIX to nf-core Modules

## Overview

This document describes the migration from local `samtools_faidx` and `tabix_vcf` modules to their nf-core equivalents. This migration improves standardization, removes local container dependencies, and adds automatic version tracking.

## Migration Summary

### SAMTOOLS_FAIDX Migration

**Status**: ✅ **COMPLETED**

**Local Module**: `modules/local/samtools_faidx.nf`  
**nf-core Module**: `modules/nf-core/samtools/faidx/main.nf`

**Files Updated**:
- `workflows/preparation/prepare_data_complete_grch37.nf`
- `workflows/preparation/prepare_data_complete_grch38.nf`
- `nextflow.config`

**Files Already Using nf-core**:
- `workflows/prepare_data_grch37.nf` (already migrated)
- `workflows/prepare_data_grch38.nf` (already migrated)

### TABIX_VCF Status

**Status**: ⚠️ **NOT USED**

The local `modules/local/tabix_vcf.nf` module was imported in `prepare_data_complete_grch37.nf` but **never actually called** in the workflow. The import statement has been removed as part of this cleanup.

---

## Technical Changes

### 1. SAMTOOLS_FAIDX Migration

#### Input Signature Change

| Aspect | Local Module | nf-core Module |
|--------|-------------|----------------|
| **Input 1** | `path fasta` | `tuple val(meta), path(fasta)` |
| **Input 2** | - | `tuple val(meta2), path(fai)` |
| **Input 3** | - | `val get_sizes` |
| **Output** | `path "${fasta}.fai"`, emit: `fai` | `tuple val(meta), path("*.fai")`, emit: `fai` |
| **Container** | Local singularity image | Public biocontainers/samtools |
| **Versions** | None | `versions.yml` output |

#### Code Changes

**Before (Local Module)**:
```groovy
include { SAMTOOLS_FAIDX } from '../../modules/local/samtools_faidx'

SAMTOOLS_FAIDX(GUNZIP.out.gunzip)
reference_fai = SAMTOOLS_FAIDX.out.fai.map { meta, file -> file }
```

**After (nf-core Module)**:
```groovy
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'

// nf-core SAMTOOLS_FAIDX requires: tuple [meta, fasta], tuple [meta2, fai], val(get_sizes)
SAMTOOLS_FAIDX(
    GUNZIP.out.gunzip,
    [[id: 'null'], []],  // Empty fai channel (we're generating it)
    false                // Don't generate .sizes file
)
reference_fai = SAMTOOLS_FAIDX.out.fai.map { meta, file -> file }
```

#### Key Parameters

**Required Inputs**:
1. **tuple [meta, fasta]**: Input FASTA file with metadata map
2. **tuple [meta2, fai]**: Existing FAI file (use empty `[[id: 'null'], []]` when generating new)
3. **val get_sizes**: Boolean flag to generate `.sizes` file (use `false` for standard FAI only)

**Outputs**:
- `tuple val(meta), path("*.fai")` - FASTA index file
- `tuple val(meta), path("*.sizes")` - Optional chromosome sizes file (if `get_sizes = true`)
- `path "versions.yml"` - Version tracking

#### Configuration Added

```groovy
withName: 'SAMTOOLS_FAIDX' {
    cpus = 1
    memory = 2.GB
    time = 1.h
    publishDir = [
        path: { "${params.references_dir}" },
        mode: 'copy',
        enabled: params.publish_faidx
    ]
}
```

**New Parameter**:
```groovy
publish_faidx = true  // Control SAMTOOLS_FAIDX output publishing
```

---

## Benefits of Migration

### 1. No Local Container Dependencies
- **Before**: Required local singularity image at `file:///singularity_images/samtools_latest.sif`
- **After**: Uses public biocontainers (`biocontainers/samtools:1.22.1--h96c455f_0`)

### 2. Automatic Version Tracking
- nf-core module generates `versions.yml` for reproducibility
- Captures exact samtools version used

### 3. Standardization
- Follows nf-core conventions with meta maps
- Consistent with other nf-core modules in pipeline
- Easier for others to understand and contribute

### 4. Easy Updates
```bash
# Update to latest nf-core version
nf-core modules update samtools/faidx
```

### 5. Additional Features
- Can generate `.sizes` file by setting `get_sizes = true`
- Supports compressed FASTA input (generates `.gzi` file)
- Customizable via `task.ext.args`

---

## Workflows Updated

### ✅ prepare_data_complete_grch37.nf
- Removed unused `TABIX_VCF` import
- Updated `SAMTOOLS_FAIDX` to nf-core module
- Added proper input parameters for nf-core signature

### ✅ prepare_data_complete_grch38.nf
- Updated `SAMTOOLS_FAIDX` to nf-core module
- Added proper input parameters for nf-core signature

### ✅ prepare_data_grch37.nf
- Already using nf-core `SAMTOOLS_FAIDX` (no changes needed)

### ✅ prepare_data_grch38.nf
- Already using nf-core `SAMTOOLS_FAIDX` (no changes needed)

---

## Testing

### Test GRCh37 Complete Data Preparation
```bash
nextflow run main.nf -profile complete_grch37 \
    --prepare_complete_data true \
    --skip_bam_download true \
    --skip_singularity_download true
```

### Test GRCh38 Complete Data Preparation
```bash
nextflow run main.nf -profile complete_grch38 \
    --prepare_complete_data true \
    --skip_bam_download true
```

### Expected Outputs

**SAMTOOLS_FAIDX**:
- `${params.references_dir}/<reference>.fasta.fai` - FASTA index file
- `work/<task_dir>/versions.yml` - Version tracking

---

## Files That Can Be Removed

After this PR is merged, the following files are no longer referenced and can be removed:

1. ✅ `modules/local/samtools_faidx.nf` - Replaced by nf-core module
2. ⚠️ `modules/local/tabix_vcf.nf` - Never used (consider removing or migrating if needed later)

**Note**: The local `tabix_vcf.nf` module was imported but never called in any workflow. If VCF indexing is needed in the future, consider using the nf-core `TABIX_TABIX` module instead.

---

## Migration Progress

- ✅ **bedtools** - Completed (PR #1)
- ✅ **gunzip** - Completed (PR #2)
- ✅ **samtools_faidx** - Completed (this PR)
- ⚠️ **tabix_vcf** - Unused (import removed)

### Remaining Local Modules

The following local modules remain in `modules/local/`:
- `create_target_beds.nf` - Custom logic specific to this pipeline
- `cutesv.nf` - SV caller wrapper
- `download_*.nf` - Custom download utilities
- `manta.nf` - SV caller wrapper
- `pbsv.nf` - SV caller wrapper
- `sniffles.nf` - SV caller wrapper
- `truvari.nf` - Benchmarking wrapper

These modules contain pipeline-specific logic and are appropriate to remain as local modules.

---

## nf-core Module Reference

### SAMTOOLS_FAIDX Documentation

**Module**: `samtools/faidx`  
**Container**: `biocontainers/samtools:1.22.1--h96c455f_0`  
**Tool Version**: samtools 1.22.1

**Command**:
```bash
samtools faidx <fasta> [args]
```

**Common Arguments** (via `task.ext.args`):
- `-o <file>` - Output to specific file
- `-f` - Force overwrite existing index

**Example Usage**:
```groovy
// Generate FAI only
SAMTOOLS_FAIDX(
    ch_fasta,              // tuple [meta, fasta]
    [[id: 'null'], []],    // Empty FAI (generating new)
    false                  // Don't generate .sizes
)

// Generate both FAI and .sizes
SAMTOOLS_FAIDX(
    ch_fasta,              // tuple [meta, fasta]
    [[id: 'null'], []],    // Empty FAI (generating new)
    true                   // Generate .sizes file
)
```

---

## Troubleshooting

### Issue: "No such file or directory" for FAI input

**Cause**: Trying to provide an existing FAI file but path is incorrect.

**Solution**: When generating a new FAI file, use empty channel:
```groovy
[[id: 'null'], []]  // Empty fai channel
```

### Issue: Output not published to references_dir

**Cause**: `params.publish_faidx` may be set to `false`.

**Solution**: Check configuration:
```groovy
params.publish_faidx = true
```

### Issue: Container not found

**Cause**: Singularity trying to pull container but no internet access.

**Solution**: Pre-pull container or use local cache:
```bash
# Pre-pull container
singularity pull docker://biocontainers/samtools:1.22.1--h96c455f_0
```

---

## Related Documentation

- [bedtools Migration](BEDTOOLS_MIGRATION.md)
- [gunzip Migration](GUNZIP_MIGRATION.md)
- [nf-core modules documentation](https://nf-co.re/modules)
- [SAMTOOLS_FAIDX nf-core module](https://github.com/nf-core/modules/tree/master/modules/nf-core/samtools/faidx)

---

## Questions?

If you have questions about this migration or encounter issues:
1. Check the nf-core module documentation
2. Review the example workflows in this repository
3. Open an issue on GitHub

**Migration completed**: 2025-11-20  
**Author**: Seqera AI
