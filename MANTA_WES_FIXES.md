# Manta WES Configuration Fixes & Improvements

## Summary
Fixed and improved the MANTA_WES process to properly support Whole Exome Sequencing analysis:

### Core Fixes
1. Added `--exome` flag to turn off depth filters for WES data
2. Added `--callRegions` parameter support to restrict variant calling to target regions

### Improvements
3. Auto-detection of tabix index files (.tbi)
4. Cleaned up container configuration (removed redundant overrides)

## Changes Made

### 1. Added `wes_sequencing_targets` parameter (nextflow.config)
```groovy
wes_sequencing_targets = null  // BED file for WES target regions (used by Manta --callRegions)
```

This parameter should point to a bgzip-compressed and tabix-indexed BED file containing the exome target regions.

### 2. Updated MANTA_WES process call (workflows/sv_calling.nf)
The process now properly passes the target BED file and its index to Manta:

```groovy
// Prepare target BED if provided
def target_bed = params.wes_sequencing_targets ? 
    file(params.wes_sequencing_targets, checkIfExists: true) : []
def target_bed_tbi = params.wes_sequencing_targets ? 
    file("${params.wes_sequencing_targets}.tbi", checkIfExists: true) : []

ch_illumina_wes_bam = Channel.value([
    [id: 'Illumina_WES', technology: 'Illumina_WES', tool: 'Manta'],
    file(params.illumina_wes_bam, checkIfExists: !is_remote),
    file("${params.illumina_wes_bam}.bai", checkIfExists: !is_remote),
    target_bed,      // Now properly set
    target_bed_tbi   // Now properly set
])
```

### 3. Added `--exome` flag to MANTA_WES configuration (conf/modules.config)
```groovy
withName: 'MANTA_WES' {
    container = 'quay.io/biocontainers/manta:1.6.0--h9ee0642_1'
    ext.args = '--exome'  // Enable WES mode: turn off depth filters
}

withName: 'MANTA_WGS' {
    container = 'quay.io/biocontainers/manta:1.6.0--h9ee0642_1'
}
```

## Improvement Details

### Auto-detection of Tabix Index (NEW!)
The pipeline now automatically detects the `.tbi` index file for the target BED:
- **No longer requires explicit index specification**
- Provides helpful warning if index is missing
- Includes instructions for creating the index

If the index is missing, you'll see:
```
⚠️  Tabix index not found for WES target regions!
    Expected: /path/to/targets.bed.gz.tbi
    
    To create the index, run:
    tabix -p bed /path/to/targets.bed.gz
```

### Container Configuration
Container specifications are explicitly defined in `conf/modules.config`:

```groovy
withName: 'MANTA_WES' {
    container = 'quay.io/biocontainers/manta:1.6.0--h9ee0642_1'
    ext.args = '--exome'  // Turn off depth filters for WES data
}
```

**Why explicit container paths?**
Some nf-core modules have incorrect or incomplete container paths (e.g., missing `quay.io/` registry prefix). The explicit specifications ensure:
- ✅ Correct registry paths (`quay.io/biocontainers`)
- ✅ Pinned versions for reproducibility
- ✅ Compatibility with both Docker and Singularity
- ✅ No runtime errors from missing containers

## Usage

To use the WES target regions, provide the path to your bgzip-compressed BED file.
The tabix index (.tbi) will be auto-detected:

### Command line:
```bash
nextflow run main.nf \
  --illumina_wes_bam /path/to/wes.bam \
  --wes_sequencing_targets /path/to/targets.bed.gz \
  --fasta /path/to/reference.fasta \
  -profile singularity
```

### In params.yaml:
```yaml
illumina_wes_bam: '/path/to/wes.bam'
wes_sequencing_targets: '/path/to/targets.bed.gz'
```

## Notes

1. **BED file format**: The target BED file MUST be:
   - bgzip-compressed (`.bed.gz`)
   - tabix-indexed (`.bed.gz.tbi` file should exist - will auto-detect)
   - Containing the same chromosome names as the reference genome

2. **Preparing the BED file**:
   ```bash
   # If you have an uncompressed BED file:
   bgzip targets.bed
   tabix -p bed targets.bed.gz
   
   # The pipeline will automatically find targets.bed.gz.tbi
   ```

3. **What --exome does**: 
   - Turns off depth-based filters that assume whole genome coverage
   - Essential for WES data to prevent false negatives in exome regions

4. **What --callRegions does**:
   - Restricts variant calling to specified target regions
   - Reduces runtime and focuses on relevant genomic regions
   - Note: The full genome is still used to estimate statistics (like expected depth)

## References

- [Manta User Guide - Exome/Targeted Sequencing](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#exometargeted)
- Manta documentation states: "For exome/targeted sequencing, use `--exome` to turn off depth filters"
