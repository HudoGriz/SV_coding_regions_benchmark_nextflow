# Migration to nf-core Modules: BGZIP and TABIX

## Summary

Successfully migrated from local `BGZIP_TABIX` and `TABIX_VCF` modules to the standardized nf-core `TABIX_BGZIPTABIX` module.

## Changes Made

### 1. Module Installation

Installed the nf-core `tabix/bgziptabix` module:
```bash
nf-core modules install tabix/bgziptabix
```

The module was added to `modules.json` with git SHA `61846a0d757ff9c0682ba3e28ab0441afd95ad7e`.

### 2. Main Workflow (`main.nf`)

**Updated imports:**
```groovy
# BEFORE:
include { BGZIP_TABIX as BGZIP_TABIX_PBSV } from './modules/local/bgzip_tabix'
include { BGZIP_TABIX as BGZIP_TABIX_CUTESV_PACBIO } from './modules/local/bgzip_tabix'
include { BGZIP_TABIX as BGZIP_TABIX_CUTESV_ONT } from './modules/local/bgzip_tabix'

# AFTER:
include { TABIX_BGZIPTABIX as BGZIP_TABIX_PBSV } from './modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as BGZIP_TABIX_CUTESV_PACBIO } from './modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIPTABIX as BGZIP_TABIX_CUTESV_ONT } from './modules/nf-core/tabix/bgziptabix/main'
```

**Updated output channel references:**
```groovy
# BEFORE:
BGZIP_TABIX_CUTESV_PACBIO.out.vcf
BGZIP_TABIX_PBSV.out.vcf
BGZIP_TABIX_CUTESV_ONT.out.vcf

# AFTER:
BGZIP_TABIX_CUTESV_PACBIO.out.gz_index
BGZIP_TABIX_PBSV.out.gz_index
BGZIP_TABIX_CUTESV_ONT.out.gz_index
```

### 3. Configuration Files

#### `conf/modules.config`

**BEFORE:**
```groovy
// BGZIP and TABIX for VCF compression and indexing
withName: 'BGZIP_TABIX' {
    container = 'quay.io/biocontainers/tabix:1.11--hdfd78af_0'
}
```

**AFTER:**
```groovy
// TABIX_BGZIPTABIX for VCF compression and indexing (nf-core module)
// Container is defined in the nf-core module
withName: 'BGZIP_TABIX.*' {
    ext.args2 = '-p vcf'
}
```

**Note:** The configuration was simplified to use a single wildcard pattern that matches all BGZIP_TABIX processes.

**Note:** Container definitions are now handled by the nf-core module itself. The nf-core module uses:
- `community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa` (Docker)
- `https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/92/92859404d861ae01afb87e2b789aebc71c0ab546397af890c7df74e4ee22c8dd/data` (Singularity)

### 4. No Changes Required

- **`nextflow.config`**: The wildcard pattern `withName: 'BGZIP_TABIX.*'` continues to work with the new process names
- **Workflow logic**: Process invocations remain the same since we used aliases
- **Output structure**: Both emit `tuple [meta, vcf.gz, vcf.gz.tbi]`, just with different emit names

## Key Differences Between Local and nf-core Modules

### Local `BGZIP_TABIX` Module
- **Output emit name**: `vcf`
- **Output structure**: `tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi")`
- **Prefix**: `${meta.id}_${meta.tool}.vcf.gz`
- **Container**: `quay.io/biocontainers/tabix:1.11--hdfd78af_0`

### nf-core `TABIX_BGZIPTABIX` Module
- **Output emit name**: `gz_index`
- **Output structure**: `tuple val(meta), path("*.gz"), path("*.{tbi,csi}")`
- **Prefix**: `${prefix}.${input.getExtension()}.gz` (preserves original extension)
- **Container**: `community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa`
- **Additional features**:
  - Multi-threading support via `task.cpus`
  - Support for both `.tbi` and `.csi` index formats
  - Configurable via `ext.args` (bgzip) and `ext.args2` (tabix)

## Benefits of Migration

1. **Standardization**: Using community-maintained modules ensures consistency across pipelines
2. **Updates**: Automatic access to bug fixes and improvements via `nf-core modules update`
3. **Multi-threading**: nf-core module supports parallel compression/indexing
4. **Flexibility**: Supports both TBI and CSI index formats
5. **Modern containers**: Uses newer htslib version (1.21 vs 1.11)

## Testing

The configuration was validated successfully:
```bash
nextflow config -profile test_nfcore
# Successfully loaded without errors
```

## Local Modules Status

The following local modules are **no longer used** and can be optionally removed:
- `modules/local/bgzip_tabix.nf`
- `modules/local/tabix_vcf.nf` (was not used anywhere)

**Note:** Keeping them for reference is fine, but they are not included in the workflow anymore.

## Verification Steps

To verify the migration works correctly:

1. Run a test with the nf-core test profile:
   ```bash
   nextflow run main.nf -profile test_nfcore,docker
   ```

2. Check that VCF files are properly compressed and indexed in the output directory:
   ```bash
   ls -lh test_results_nfcore/calls/*/sv/*/
   ```

3. Verify the index files can be used with tabix:
   ```bash
   tabix -l test_results_nfcore/calls/*/sv/*/*.vcf.gz
   ```

## Migration Date

Date: 2025-11-25
Nextflow version: 25.04.7
nf-core/tools version: 3.3.2

## Additional Notes

- The `conf/modules_docker.config` file was removed in favor of a simplified configuration
- All container definitions are now managed by the nf-core module itself
- The wildcard pattern `BGZIP_TABIX.*` in `nextflow.config` continues to work seamlessly
