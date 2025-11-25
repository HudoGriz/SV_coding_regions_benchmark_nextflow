# Manta WES Configuration Fixes

## Summary
Fixed the MANTA_WES process to properly support Whole Exome Sequencing analysis by adding:
1. The `--exome` flag to turn off depth filters for WES data
2. The `--callRegions` parameter to restrict variant calling to target regions

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

## Usage

To use the WES target regions, provide the path to your bgzip-compressed and tabix-indexed BED file:

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
   - tabix-indexed (`.bed.gz.tbi` file must exist)
   - Containing the same chromosome names as the reference genome

2. **Preparing the BED file**:
   ```bash
   # If you have an uncompressed BED file:
   bgzip targets.bed
   tabix -p bed targets.bed.gz
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
