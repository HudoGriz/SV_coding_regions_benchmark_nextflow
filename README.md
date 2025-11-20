# SV Calling and Benchmarking Pipeline

A Nextflow DSL2 pipeline for calling structural variants across multiple sequencing technologies and benchmarking results with Truvari.

## Overview

This pipeline processes BAM files from different sequencing technologies to call structural variants and benchmark the calls against a truth set:

**Supported Technologies & Tools:**
- **Illumina WES**: Manta
- **Illumina WGS**: Manta
- **PacBio**: CuteSV, Pbsv
- **ONT**: CuteSV, Sniffles

**Benchmarking:**
- Truvari benchmarking against truth set VCF
- Three target sets: high confidence, gene panel, WES UTR regions

## Requirements

- Nextflow >= 23.04.0
- Singularity/Apptainer
- Local Singularity images in `singularity_images/` directory:
  - `manta_latest.sif`
  - `cutesv_latest.sif`
  - `pbsv_latest.sif`
  - `sniffles_latest.sif`
  - `samtools_latest.sif`
  - `truvari_modded.sif`

## Quick Start

### 1. Prepare Input Data

Ensure you have:
- BAM files (indexed with .bai)
- Reference genome FASTA (indexed with .fai)
- Benchmarking VCF (indexed with .tbi)
- Target BED files

### 2. Configure Parameters

Edit `params.yaml` with your file paths:

```yaml
fasta: '/path/to/reference.fasta'
illumina_wes_bam: '/path/to/illumina_wes.bam'
benchmark_vcf: '/path/to/truth_set.vcf.gz'
# ... etc
```

### 3. Run the Pipeline

```bash
# Run with all technologies
nextflow run main.nf -params-file params.yaml -profile singularity

# Run specific technologies only (comment out unwanted BAMs in params.yaml)
nextflow run main.nf -params-file params.yaml -profile singularity

# Resume a previous run
nextflow run main.nf -params-file params.yaml -profile singularity -resume
```

## Pipeline Structure

```
sv-calling-pipeline/
├── main.nf                 # Main workflow file
├── nextflow.config         # Configuration file
├── params.yaml            # Example parameters
├── modules/
│   └── local/
│       ├── manta.nf       # Illumina SV calling (Manta)
│       ├── cutesv.nf      # Long-read SV calling (CuteSV)
│       ├── pbsv.nf        # PacBio SV calling (Pbsv)
│       ├── sniffles.nf    # ONT SV calling (Sniffles)
│       └── truvari.nf     # Benchmarking with Truvari
└── singularity_images/    # Local Singularity containers
```

## Output Structure

```
results/
├── calls/
│   ├── Illumina_WES/
│   │   └── sv/manta/
│   ├── Illumina_WGS/
│   │   └── sv/manta/
│   ├── PacBio/
│   │   ├── sv/cutesv/
│   │   └── sv/pbsv/
│   └── ONT/
│       ├── sv/cutesv/
│       └── sv/sniffles/
├── real_intervals/
│   ├── Illumina_WES/
│   │   └── truvari_benchmark/
│   ├── Illumina_WGS/
│   │   └── truvari_benchmark/
│   ├── PacBio/
│   │   └── truvari_benchmark/
│   └── ONT/
│       └── truvari_benchmark/
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    ├── execution_trace.txt
    └── pipeline_dag.svg
```

## Parameters

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `fasta` | Reference genome FASTA file |
| `benchmark_vcf` | Truth set VCF for benchmarking |
| `*_targets` | BED files defining target regions |

### Input BAM Files (Optional)

| Parameter | Description |
|-----------|-------------|
| `illumina_wes_bam` | Illumina WES BAM file |
| `illumina_wgs_bam` | Illumina WGS BAM file |
| `pacbio_bam` | PacBio BAM file |
| `ont_bam` | ONT BAM file |

**Note:** Only provide BAM files for technologies you want to analyze. The pipeline will automatically skip technologies without BAM files.

### Truvari Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `truvari_refdist` | 500 | Reference distance for matching |
| `truvari_pctsize` | 0.7 | Percent size similarity |
| `truvari_pctovl` | 0 | Percent overlap |
| `truvari_pctseq` | 0 | Percent sequence similarity |

WES-specific parameters can be set with `truvari_wes_*` prefix.

### Resource Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_cpus` | 24 | Maximum CPUs per process |
| `max_memory` | 128.GB | Maximum memory per process |
| `max_time` | 48.h | Maximum time per process |

## Profiles

- `singularity`: Use Singularity containers (recommended)
- `local`: Run on local machine
- `slurm`: Submit jobs to SLURM cluster
- `standard`: Default profile (Singularity + local executor)

## Advanced Usage

### Running Only Specific Technologies

Comment out unwanted BAM files in `params.yaml`:

```yaml
# illumina_wes_bam: '/path/to/wes.bam'  # Commented out - will skip
illumina_wgs_bam: '/path/to/wgs.bam'    # Will run
```

### Overriding Parameters on Command Line

```bash
nextflow run main.nf \
  -params-file params.yaml \
  --max_cpus 48 \
  --outdir custom_output
```

### Using Custom Container Paths

Edit `nextflow.config`:

```groovy
params {
    container_manta = "file:///custom/path/manta.sif"
    // ... etc
}
```

## Key Differences from Bash Pipeline

1. **Parallel Execution**: All technologies run in parallel automatically
2. **Resume Capability**: Use `-resume` to continue from last completed step
3. **Resource Management**: Automatic CPU/memory allocation per process
4. **Scalability**: Easy to run on HPC clusters with profile changes
5. **Reproducibility**: Explicit container and parameter versioning
6. **No Manual Checkpointing**: Nextflow handles execution tracking

## Troubleshooting

### Singularity Mount Issues

If containers can't access files, ensure paths are accessible:
- Nextflow automatically mounts paths (v23.09+)
- Use absolute paths in parameters
- Check Singularity configuration: `singularity.autoMounts = true`

### BAM Index Not Found

Ensure `.bai` files exist alongside BAM files:
```bash
samtools index your_file.bam
```

### Out of Memory Errors

Increase memory in `nextflow.config`:
```groovy
withName: 'PROCESS_NAME' {
    memory = { check_max(64.GB * task.attempt, 'memory') }
}
```

## Citation

If you use this pipeline, please cite:
- Nextflow: https://doi.org/10.1038/nbt.3820
- Individual tools (Manta, CuteSV, Pbsv, Sniffles, Truvari)

## License

This pipeline is provided as-is for research purposes.
