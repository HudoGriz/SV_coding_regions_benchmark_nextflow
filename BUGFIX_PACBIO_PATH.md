# 🐛 Bug Fix: PacBio Test Data File Path

## Issue

The test profile was failing with a 404 error when trying to download PacBio test data:

```
ERROR ~ Error executing process > 'CUTESV_PACBIO'
Caused by:
  Can't stage file https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/pacbio/bam/test_hifi.sorted.bam
  -- reason: Unable to access path
```

## Root Cause

The configuration was referencing a **non-existent file**: `test_hifi.sorted.bam`

**File that doesn't exist:**
- ❌ `data/genomics/homo_sapiens/pacbio/bam/test_hifi.sorted.bam`

**File that actually exists:**
- ✅ `data/genomics/homo_sapiens/pacbio/bam/test.sorted.bam`

## Investigation

Checked the nf-core/test-datasets repository (modules branch):

```bash
# Available PacBio BAM files in nf-core/test-datasets
curl -s "https://api.github.com/repos/nf-core/test-datasets/contents/data/genomics/homo_sapiens/pacbio/bam?ref=modules"
```

**Available files:**
- test.sorted.bam ✅
- test.sorted.bam.bai ✅
- alz.bam
- alz.ccs.bam
- NA03697B2_downsampled.pbmm2.repeats.bam
- NA037562_downsampled.pbmm2.repeats.bam
- (and others...)

**Missing file:**
- test_hifi.sorted.bam ❌

## Fix Applied

### 1. Updated Configuration File

**File:** `conf/test_nfcore.config`

```diff
- pacbio_bam = "${test_data_base}/genomics/homo_sapiens/pacbio/bam/test_hifi.sorted.bam"
+ pacbio_bam = "${test_data_base}/genomics/homo_sapiens/pacbio/bam/test.sorted.bam"
```

### 2. Updated Documentation

**File:** `docs/TESTING_LONG_READ_CALLERS.md`

```diff
- **PacBio HiFi**: `data/genomics/homo_sapiens/pacbio/bam/test_hifi.sorted.bam`
+ **PacBio**: `data/genomics/homo_sapiens/pacbio/bam/test.sorted.bam`
```

Also corrected Illumina WGS filename from `test.paired_end.sorted.bam` to `test2.paired_end.sorted.bam`.

## Verification

All test data files now verified to exist:

```bash
# Illumina WES
curl -I https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test.paired_end.sorted.bam
# ✅ HTTP/2 200

# Illumina WGS
curl -I https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/illumina/bam/test2.paired_end.sorted.bam
# ✅ HTTP/2 200

# PacBio BAM
curl -I https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/pacbio/bam/test.sorted.bam
# ✅ HTTP/2 200

# PacBio BAM index
curl -I https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/pacbio/bam/test.sorted.bam.bai
# ✅ HTTP/2 200

# ONT BAM
curl -I https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/nanopore/bam/test.sorted.bam
# ✅ HTTP/2 200

# Reference genome
curl -I https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/genome/genome.fasta
# ✅ HTTP/2 200
```

## Commit Details

**Commit Hash:** `192c3a7`  
**Branch:** `seqera-ai/20251120-132138-add-ci-testing`  
**PR:** #7

**Files Changed:**
- `conf/test_nfcore.config` (fixed PacBio path)
- `docs/TESTING_LONG_READ_CALLERS.md` (updated documentation)

## Expected Behavior After Fix

The test profile should now successfully:
1. ✅ Download PacBio test BAM file
2. ✅ Run CUTESV_PACBIO process
3. ✅ Run PBSV_DISCOVER process
4. ✅ Run PBSV_CALL process
5. ✅ Complete all PacBio SV calling steps

## Testing

To verify the fix locally:

```bash
# Test just the PacBio data download
nextflow run . -profile test_nfcore,docker --outdir test_results

# Check that PacBio processes execute
# Should see:
# [xx/xxxxxx] CUTESV_PACBIO (PacBio)      | 1 of 1 ✔
# [xx/xxxxxx] PBSV_DISCOVER (PacBio)      | 1 of 1 ✔
# [xx/xxxxxx] PBSV_CALL (PacBio)          | 1 of 1 ✔
```

## Related Issues

This was likely a naming assumption based on the fact that the data is PacBio HiFi data, but the actual filename in the repository is just `test.sorted.bam` without the `_hifi` suffix.

## Status

✅ **Fixed and pushed to PR #7**  
🤖 **CI will test automatically on next push**  
📝 **Documentation updated**

---

**The pipeline should now successfully test all SV callers including PacBio!** 🎉
