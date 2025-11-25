# Pipeline Cleanup Audit Report
**Date:** 2025-11-25  
**Purpose:** Comprehensive audit to identify and remove unused components

## Executive Summary

### ❌ UNUSED COMPONENTS IDENTIFIED

#### 1. **Unused Local Modules** (2 files)
- `modules/local/bgzip_tabix.nf` - NEVER included/used (use nf-core TABIX_BGZIPTABIX instead)
- `modules/local/tabix_vcf.nf` - NEVER included/used

#### 2. **Unused Workflow Files** (2 files)
- `workflows/prepare_data_grch37.nf` - OLD version, replaced by `workflows/preparation/prepare_data_complete_grch37.nf`
- `workflows/prepare_data_grch38.nf` - OLD version, replaced by `workflows/preparation/prepare_data_complete_grch38.nf`

#### 3. **Disabled Config File** (1 file)
- `conf/test_simulation.config.disabled` - File already disabled by .disabled extension

---

## Detailed Analysis

### Parameters Audit
**Parameters actually USED in main.nf (19 total):**
✅ `params.help` - Display help
✅ `params.prepare_giab_resources` - Data preparation mode
✅ `params.prepare_complete_data` - Full data download mode
✅ `params.genome` - Reference genome selection
✅ `params.project_dir` - Project directory
✅ `params.skip_singularity_download` - Skip Singularity downloads
✅ `params.skip_bam_download` - Skip BAM downloads
✅ `params.skip_reference_download` - Skip reference downloads
✅ `params.download_grch37_liftover` - GRCh37 liftover option
✅ `params.outdir` - Output directory
✅ `params.illumina_wes_bam` - Illumina WES input
✅ `params.illumina_wgs_bam` - Illumina WGS input
✅ `params.pacbio_bam` - PacBio input
✅ `params.ont_bam` - ONT input
✅ `params.skip_pbsv` - Skip PBSV tool
✅ `params.benchmark_vcf` - Truth set VCF
✅ `params.skip_benchmarking` - Skip benchmarking
✅ `params.simulate_targets` - Enable simulation
✅ `params.num_simulations` - Number of simulations
✅ `params.gencode_gtf` - GTF annotation
✅ `params.gather_statistics` - Generate statistics

**Parameters in nextflow_schema.json NOT used in code:**
All schema parameters are now properly aligned with actual usage after our recent schema update.

---

### Module Files Audit

#### ✅ **USED Local Modules (9 files):**
1. `modules/local/create_target_beds.nf` - Used in preparation workflows
2. `modules/local/download_annotations.nf` - Used in preparation workflows
3. `modules/local/download_bam.nf` - Used in preparation workflows
4. `modules/local/download_reference.nf` - Used in preparation workflows
5. `modules/local/download_singularity.nf` - Used in GRCh37 preparation
6. `modules/local/download_truth_set.nf` - Used in preparation workflows
7. `modules/local/gather_statistics.nf` - Used in analysis_and_plots workflow
8. `modules/local/simulate_targets.nf` - Used in simulate_and_benchmark workflow

#### ❌ **UNUSED Local Modules (2 files):**
1. `modules/local/bgzip_tabix.nf` - **NEVER INCLUDED** anywhere
   - Replacement: Use nf-core `TABIX_BGZIPTABIX` module instead
   - Already using: `TABIX_BGZIPTABIX as BGZIP_TABIX_PBSV` in sv_calling.nf

2. `modules/local/tabix_vcf.nf` - **NEVER INCLUDED** anywhere
   - Likely leftover from old code
   - Similar functionality available in nf-core modules

---

### Workflow Files Audit

#### ✅ **USED Workflows (8 files):**
1. `workflows/prepare_references.nf` - ✅ Included in main.nf
2. `workflows/sv_calling.nf` - ✅ Included in main.nf
3. `workflows/benchmarking.nf` - ✅ Included in main.nf
4. `workflows/simulate_and_benchmark.nf` - ✅ Included in main.nf
5. `workflows/analysis_and_plots.nf` - ✅ Included in main.nf
6. `workflows/prepare_giab_resources.nf` - ✅ Included in main.nf
7. `workflows/preparation/prepare_data_complete_grch37.nf` - ✅ Included in main.nf
8. `workflows/preparation/prepare_data_complete_grch38.nf` - ✅ Included in main.nf

#### ❌ **UNUSED Workflows (2 files):**
1. `workflows/prepare_data_grch37.nf` - **NEVER INCLUDED** in main.nf
   - OLD VERSION
   - Replaced by: `workflows/preparation/prepare_data_complete_grch37.nf`
   - Contains outdated module includes

2. `workflows/prepare_data_grch38.nf` - **NEVER INCLUDED** in main.nf
   - OLD VERSION
   - Replaced by: `workflows/preparation/prepare_data_complete_grch38.nf`
   - Contains outdated module includes

---

### Configuration Files Audit

#### ✅ **USED Configs (5 files):**
1. `nextflow.config` - Main config (includes all others)
2. `conf/base.config` - Base process resources
3. `conf/modules.config` - Module-specific containers
4. `conf/modules_docker.config` - Docker-specific settings
5. `conf/test.config` - Test profile
6. `conf/test_nfcore.config` - nf-core test profile

#### ⚠️ **DISABLED Config (1 file):**
1. `conf/test_simulation.config.disabled` - Already disabled by extension
   - Not included in nextflow.config
   - Can be safely removed

---

## Cleanup Actions Recommended

### 🗑️ **Files to DELETE (5 total):**

1. **Unused Local Modules (2):**
   ```bash
   rm modules/local/bgzip_tabix.nf
   rm modules/local/tabix_vcf.nf
   ```

2. **Obsolete Workflows (2):**
   ```bash
   rm workflows/prepare_data_grch37.nf
   rm workflows/prepare_data_grch38.nf
   ```

3. **Disabled Config (1):**
   ```bash
   rm conf/test_simulation.config.disabled
   ```

---

## Verification Steps After Cleanup

1. ✅ Check that all includes in workflow files resolve correctly
2. ✅ Verify no references to deleted files exist
3. ✅ Run basic pipeline validation
4. ✅ Ensure documentation doesn't reference removed files

---

## Benefits of Cleanup

- **Reduced confusion** - No old/duplicate files
- **Clearer structure** - Only active code present
- **Easier maintenance** - Less to review and update
- **Smaller repository** - Faster clones and checkouts

---

## Notes

- All nf-core modules are properly used and should be retained
- The `preparation/` subdirectory structure is clean and well-organized
- Module container definitions in `conf/modules.config` match actual processes
- No orphaned process definitions found in configs
