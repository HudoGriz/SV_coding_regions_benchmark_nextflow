# Bash to Nextflow Conversion Notes

## Overview

This document describes how the original bash scripts were converted to Nextflow DSL2.

## Original Bash Pipeline Structure

The original pipeline consisted of several bash scripts:

1. **analysis.sh** - Main orchestrator script
2. **Illumina_wes.sh** - Illumina WES SV calling
3. **Illumina_wgs.sh** - Illumina WGS SV calling  
4. **Pacbio_wgs.sh** - PacBio SV calling
5. **ONT_wgs.sh** - ONT SV calling
6. **seq_technologies.sh** - Benchmarking coordinator
7. **target_sets.sh** - Individual benchmarking runs
8. **config.sh** - Configuration variables
9. **params_truvari.sh** - Truvari parameters

## Conversion Mapping

### Bash Script → Nextflow Component

| Bash Script | Nextflow Component | Notes |
|-------------|-------------------|-------|
| `analysis.sh` | `main.nf` workflow block | Main orchestration |
| `Illumina_*.sh` | `modules/local/manta.nf` | Process for Manta |
| `Pacbio_wgs.sh` | `modules/local/cutesv.nf` + `pbsv.nf` | Split into two processes |
| `ONT_wgs.sh` | `modules/local/cutesv.nf` + `sniffles.nf` | CuteSV reused |
| `target_sets.sh` | `modules/local/truvari.nf` | Single process, multiple invocations |
| `config.sh` | `nextflow.config` + `params.yaml` | Parameters and configuration |
| `params_truvari.sh` | `params` section in config | Truvari parameters |

### Manual Checkpointing → Nextflow Resumability

**Bash Approach:**
```bash
LOGFILE="$project_dir/$run_name/pipeline_status.log"
is_done() { grep -qx "$1" "$LOGFILE" 2>/dev/null }
mark_done() { echo "$1" >> "$LOGFILE" }
```

**Nextflow Approach:**
```bash
# Just use -resume flag
nextflow run main.nf -resume
```

Nextflow automatically tracks completed tasks and resumes from failures.

### Sequential Execution → Parallel DAG

**Bash Approach:**
- Sequential steps with manual `is_done` checks
- Parallel execution via `&` and `wait`
- Manual dependency management

**Nextflow Approach:**
- Automatic parallel execution
- Implicit DAG based on channel dependencies
- All SV calling runs in parallel automatically
- Benchmarking waits for all SV calling to complete

### Container Binding

**Bash Approach:**
```bash
singularity exec \
    -B $project_dir \
    -B $analysis_dir/sv/manta:/out \
    $project_dir/singularity_images/manta_latest.sif \
    command
```

**Nextflow Approach:**
```groovy
// In nextflow.config
singularity.autoMounts = true
process.container = "file:///path/to/image.sif"

// Nextflow handles all mounting automatically
```

### Variable Substitution

**Bash Approach:**
```bash
source $project_dir/scripts/configs/config.sh
export threads=24
command --threads $threads
```

**Nextflow Approach:**
```groovy
// In process
script:
"""
command --threads ${task.cpus}
"""
```

## Key Improvements in Nextflow Version

### 1. **True Parallelization**

**Bash**: Limited parallelization with `&` and `wait`
```bash
bash script1.sh &
bash script2.sh &
wait
```

**Nextflow**: Automatic parallel execution across all available resources
```groovy
// Automatically runs in parallel if resources available
MANTA_WES(...)
MANTA_WGS(...)
CUTESV_PACBIO(...)
```

### 2. **Automatic Resume**

**Bash**: Manual checkpoint system required
**Nextflow**: Built-in with `-resume` flag

### 3. **Resource Management**

**Bash**: Fixed resource allocation
```bash
export threads=24
```

**Nextflow**: Dynamic resource allocation with limits
```groovy
process {
    cpus = { check_max(8, 'cpus') }
    memory = { check_max(32.GB * task.attempt, 'memory') }
}
```

### 4. **Scalability**

**Bash**: Local execution only (or complex job submission wrapper)
**Nextflow**: Simple profile switching
```bash
# Local
nextflow run main.nf -profile local

# SLURM cluster
nextflow run main.nf -profile slurm

# AWS Batch
nextflow run main.nf -profile awsbatch
```

### 5. **Error Handling**

**Bash**: `set -e` exits on first error
**Nextflow**: Configurable retry strategies
```groovy
process {
    errorStrategy = 'retry'
    maxRetries = 2
}
```

### 6. **Provenance & Reporting**

**Bash**: No built-in tracking
**Nextflow**: Automatic generation of:
- Execution timeline
- Resource usage reports  
- Pipeline DAG visualization
- Execution trace

## Channel Design Patterns Used

### 1. **Value Channels for Single Inputs**
```groovy
ch_fasta = Channel.value(file(params.fasta))
```

### 2. **Mix for Combining Results**
```groovy
ch_all_vcfs = Channel.empty()
    .mix(MANTA_WES.out.vcf)
    .mix(MANTA_WGS.out.vcf)
```

### 3. **Combine for Cartesian Products**
```groovy
ch_benchmark_input = ch_all_vcfs
    .combine(ch_targets)
```

### 4. **Conditional Channel Creation**
```groovy
if (params.illumina_wes_bam) {
    ch_illumina_wes_bam = Channel.value(...)
}
```

## Parameter Handling

### Bash Configuration Files

```bash
# config.sh
export fasta=$project_dir/data/references/human_hs37d5.fasta
export threads=24

# In script
source $project_dir/scripts/configs/config.sh
```

### Nextflow Parameters

```groovy
// nextflow.config
params {
    fasta = null
    max_cpus = 24
}

// params.yaml
fasta: '/path/to/reference.fasta'
max_cpus: 48
```

**Benefits:**
- Type checking
- Default values
- Easy overrides
- Validation possible

## Output Organization

### Bash Approach
```bash
analysis_dir=$1
mkdir -p $analysis_dir/sv/manta
```

### Nextflow Approach
```groovy
publishDir "${params.outdir}/calls/${meta.technology}/sv/manta", 
    mode: 'copy'
```

**Benefits:**
- Automatic directory creation
- Metadata-based organization
- Multiple publish strategies (copy, symlink, move)
- Conditional publishing

## Container Management

### Bash Approach
- Hardcoded paths
- Manual binding
- Each script specifies container

### Nextflow Approach
```groovy
// Centralized in config
params {
    container_manta = "file:///path/to/manta.sif"
}

process {
    withName: 'MANTA.*' {
        container = params.container_manta
    }
}
```

**Benefits:**
- Single source of truth
- Easy to update versions
- Process-specific containers
- Container caching

## Testing Strategy

### Bash Pipeline
- Manual testing
- Full data required
- Hard to test individual steps

### Nextflow Pipeline
- Use `-stub-run` for dry runs
- Test with small datasets
- Test individual processes
- Use `-resume` for iterative development

```bash
# Dry run to test workflow logic
nextflow run main.nf -stub-run

# Test with small data
nextflow run main.nf --illumina_wes_bam small_test.bam -resume
```

## Migration Checklist

For users migrating from bash to Nextflow:

- [ ] Install Nextflow (>= 23.04.0)
- [ ] Ensure Singularity images are accessible
- [ ] Create `params.yaml` with your paths
- [ ] Verify BAM files have `.bai` indexes
- [ ] Verify FASTA has `.fai` index
- [ ] Verify VCF has `.tbi` index
- [ ] Test with single technology first
- [ ] Enable all technologies once validated
- [ ] Compare results with bash pipeline
- [ ] Set up for cluster execution if needed

## Performance Considerations

### Bash Pipeline
- Manual parallelization limited
- Sequential by default
- No automatic load balancing

### Nextflow Pipeline
- Automatic parallelization
- Efficient resource usage
- Can scale to HPC/cloud
- Automatically retries failed tasks

### Expected Speedup
- **Single machine**: Similar to bash with `&` parallelization
- **HPC cluster**: Significantly faster (10-100x depending on resources)
- **Resume capability**: Much faster for reruns after failures

## Troubleshooting

### Common Issues When Converting

1. **Path Issues**
   - Use absolute paths in `params.yaml`
   - Verify Singularity can access paths

2. **Container Differences**
   - Bash uses `-B` for binding
   - Nextflow auto-mounts (v23.09+)

3. **Environment Variables**
   - Don't rely on shell environment
   - Put everything in params

4. **File Extensions**
   - Nextflow expects companion files (.bai, .fai, .tbi)
   - Explicitly specify in channels if non-standard

## Further Optimization

Possible enhancements:

1. **Sample Sheet Input**: Replace individual BAM params with CSV samplesheet
2. **Module Imports**: Use nf-core modules where available
3. **MultiQC**: Add MultiQC report aggregation
4. **Validation**: Add parameter validation schema
5. **Stub Runs**: Add stub sections for testing
6. **Subworkflows**: Group related processes

## Conclusion

The Nextflow version provides:
- ✅ Better parallelization
- ✅ Automatic resume capability
- ✅ Easier scalability
- ✅ Better provenance tracking
- ✅ More maintainable code
- ✅ Container management
- ✅ Resource optimization

While maintaining:
- ✅ Same scientific outputs
- ✅ Same tool versions (via containers)
- ✅ Same parameter settings
