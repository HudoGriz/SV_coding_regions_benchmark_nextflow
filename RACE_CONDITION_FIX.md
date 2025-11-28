# Race Condition Fix for GATHER_STATISTICS Process

## Problem Statement

The GATHER_STATISTICS process was experiencing a race condition where it would fail with "file not found" errors when trying to read Truvari benchmark results. This occurred because:

1. R scripts read data from `params.outdir` (the published results directory)
2. TRUVARI_BENCH completes and emits output channels
3. GATHER_STATISTICS starts immediately  
4. Publishing happens asynchronously in the background
5. R scripts try to read files that haven't been published yet → **FAILURE**

## Root Cause

The pipeline was relying on Nextflow's asynchronous publishing mechanism, which has no guarantees about when files will be available in the published directory. The R scripts expected a specific directory structure:
- `real_intervals/{technology}/truvari/{tool}/{target}/`
- `simulated_intervals/benchmarks/{technology}/{tool}/{target}/`

## Solution

**Stage actual work directory files directly into GATHER_STATISTICS**, eliminating dependency on published outputs.

### Architecture Changes

#### 1. ANALYSIS_AND_PLOTS Workflow (`workflows/analysis_and_plots.nf`)

**Before:**
```groovy
workflow ANALYSIS_AND_PLOTS {
    take:
    ch_truvari_results      // [meta, summary.json]
    ch_run_dir              // String path to params.outdir

    main:
    ch_ready = ch_truvari_results
        .map { meta, summary -> 1 }
        .collect()
        .map { it.size() }

    GATHER_STATISTICS(
        ch_ready,
        ch_run_dir
    )
}
```

**After:**
```groovy
workflow ANALYSIS_AND_PLOTS {
    take:
    ch_truvari_results      // [meta, summary.json]

    main:
    ch_organized_results = ch_truvari_results
        .map { meta, summary ->
            def result_dir = summary.parent  // Get work directory
            def target_set = meta.target_set ?: 'real'
            def stage_path = target_set == 'simulated' ? 
                "simulated_intervals/benchmarks/${meta.technology}/${meta.tool}/${meta.target}" :
                "real_intervals/${meta.technology}/truvari/${meta.tool}/${meta.target}"
            
            tuple(stage_path, result_dir)
        }
        .collect()
        .map { list_of_tuples ->
            def paths = list_of_tuples.collect { it[0] }
            def dirs = list_of_tuples.collect { it[1] }
            tuple(paths, dirs)
        }

    GATHER_STATISTICS(
        ch_organized_results
    )
}
```

**Key Changes:**
- Removed `ch_run_dir` parameter
- Extract parent directory from `summary.json` (contains all Truvari outputs)
- Collect all result directories with their expected staging paths
- Pass as tuple of (paths list, directories list)

#### 2. GATHER_STATISTICS Process (`modules/local/gather_statistics.nf`)

**Before:**
```groovy
process GATHER_STATISTICS {
    input:
    val ready
    val run_dir

    script:
    """
    Rscript ${projectDir}/bin/R/paper_plots.R ${run_dir} plots tables
    Rscript ${projectDir}/bin/R/general_statistics.R ${run_dir} summary_statistics.txt
    """
}
```

**After:**
```groovy
process GATHER_STATISTICS {
    input:
    tuple val(stage_paths), path(result_dirs, stageAs: 'staged_*')

    script:
    """
    # Create directory structure
    mkdir -p real_intervals simulated_intervals/benchmarks
    
    # Stage all Truvari result directories using symlinks
    count=0
    for staged_dir in staged_*; do
        if [ -d "\$staged_dir" ]; then
            target_path=\$(echo '${stage_paths.join('\n')}' | sed -n "\$((count+1))p")
            mkdir -p "\$(dirname "\$target_path")"
            ln -sf "\$(pwd)/\$staged_dir" "\$target_path"
            count=\$((count+1))
        fi
    done
    
    # Run R scripts on staged data
    Rscript ${projectDir}/bin/R/paper_plots.R \$(pwd) plots tables
    Rscript ${projectDir}/bin/R/general_statistics.R \$(pwd) summary_statistics.txt
    """
}
```

**Key Changes:**
- Input: `tuple val(stage_paths), path(result_dirs, stageAs: 'staged_*')`
  - `stage_paths`: List of target paths (e.g., "real_intervals/ONT/truvari/cuteSV/gene_panel")
  - `result_dirs`: Actual work directories staged as `staged_0`, `staged_1`, etc.
- Script creates expected directory structure with symlinks
- R scripts read from process work directory instead of `params.outdir`

#### 3. Main Workflow Update (`main.nf`)

**Before:**
```groovy
if (params.gather_statistics && ...) {
    ch_run_dir = Channel.value(params.outdir)
    
    ANALYSIS_AND_PLOTS(
        ch_truvari_results,
        ch_run_dir
    )
}
```

**After:**
```groovy
if (params.gather_statistics && ...) {
    ANALYSIS_AND_PLOTS(
        ch_truvari_results
    )
}
```

**Key Changes:**
- Removed `ch_run_dir` channel creation
- ANALYSIS_AND_PLOTS now only needs `ch_truvari_results`

## Benefits

1. **Eliminates Race Condition**: Files are guaranteed to exist when GATHER_STATISTICS runs
2. **No Timing Dependencies**: Doesn't rely on publishing being complete
3. **Maintains Structure**: R scripts see the same directory structure they expect
4. **Efficient**: Uses symlinks instead of copying data
5. **Atomic**: All inputs are staged before process execution begins

## Testing

To verify the fix works:

```bash
# Run with benchmarking and statistics enabled
nextflow run main.nf \
    --fasta reference.fa \
    --illumina_wes_bam sample.bam \
    --benchmark_vcf truth.vcf.gz \
    --gene_panel_targets regions.bed \
    --gather_statistics \
    -profile docker
```

Check that:
1. GATHER_STATISTICS completes without "file not found" errors
2. Plots and tables are generated correctly
3. Both real and simulated intervals are processed (if applicable)

## Implementation Notes

- The `stageAs: 'staged_*'` directive ensures each directory gets a unique name
- Bash script iterates through staged directories in order
- Symlinks preserve directory structure expected by R scripts
- Process working directory serves as the "run_dir" for R scripts
