# SV Calling and Benchmarking Pipeline

A Nextflow DSL2 pipeline for calling structural variants across multiple sequencing technologies and benchmarking results with Truvari.

## Features

- **Multi-technology SV calling**: Supports Illumina (WES/WGS), PacBio, and ONT sequencing
- **Multiple callers**: Manta, CuteSV, Pbsv, Sniffles
- **Comprehensive benchmarking**: Truvari comparison against truth sets
- **Flexible target regions**: High confidence regions, gene panels, or WES with UTR regions
- **Parallel execution**: All technologies and targets processed simultaneously
- **Resume capability**: Continue from last successful step with `-resume`

## Requirements

- Nextflow >= 23.04.0
- Container engine: Singularity/Apptainer, Docker, or Conda
- Input data:
  - BAM files with indexes (.bai)
  - Reference genome FASTA with index (.fai)
  - Truth/baseline VCF with index (.tbi)
  - Target BED files (optional)

## Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/HudoGriz/SV_coding_regions_benchmark_nextflow.git
cd SV_coding_regions_benchmark_nextflow

# 2. Create a parameters file
cat > my_params.yaml << EOF
# Required inputs
fasta: /path/to/reference.fasta
vcf_baseline: /path/to/truth_set.vcf.gz

# Optional BAM files (provide only those you want to analyze)
illumina_wes_bam: /path/to/illumina_wes.bam
illumina_wgs_bam: /path/to/illumina_wgs.bam
pacbio_bam: /path/to/pacbio.bam
ont_bam: /path/to/ont.bam

# Target regions (at least one recommended)
high_confidence_targets: /path/to/high_conf.bed
gene_panel_targets: /path/to/gene_panel.bed
wes_with_utr_targets: /path/to/wes_utr.bed

# Output configuration
outdir: ./results
run_name: my_sv_analysis
EOF

# 3. Run the pipeline
nextflow run main.nf -params-file my_params.yaml -profile singularity

# 4. Resume if interrupted
nextflow run main.nf -params-file my_params.yaml -profile singularity -resume
```

## Pipeline Overview

### Supported Technologies and Tools

| Technology | SV Callers |
|-----------|-----------|
| Illumina WES | Manta |
| Illumina WGS | Manta |
| PacBio HiFi | CuteSV, Pbsv |
| ONT | CuteSV, Sniffles |

### Workflow

1. **SV Calling**: Process each BAM file with appropriate caller(s)
2. **Benchmarking**: Compare called SVs against truth set using Truvari
3. **Target Analysis**: Evaluate performance across different target regions

## Input Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `fasta` | Reference genome FASTA file |
| `vcf_baseline` | Truth/baseline VCF for benchmarking |

### BAM Files (at least one required)

| Parameter | Description |
|-----------|-------------|
| `illumina_wes_bam` | Illumina whole-exome sequencing BAM |
| `illumina_wgs_bam` | Illumina whole-genome sequencing BAM |
| `pacbio_bam` | PacBio HiFi long-read BAM |
| `ont_bam` | Oxford Nanopore BAM |

**Note**: Only provide BAM files for technologies you want to analyze.

### Target Regions (optional but recommended)

| Parameter | Description |
|-----------|-------------|
| `high_confidence_targets` | High-confidence genomic regions BED file |
| `gene_panel_targets` | Gene panel regions BED file |
| `wes_with_utr_targets` | WES capture regions including UTRs BED file |

### Optional Annotations

| Parameter | Description | Used By |
|-----------|-------------|---------|
| `tandem_repeats` | Tandem repeat regions BED | Sniffles (ONT) |
| `wes_sequencing_targets` | WES capture targets BED.gz + .tbi | Manta (WES) |

### Truvari Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `truvari_refdist` | 500 | Maximum reference distance (bp) |
| `truvari_pctsize` | 0.7 | Minimum size similarity (0-1) |
| `truvari_pctseq` | 0.0 | Minimum sequence similarity (0-1) |
| `truvari_pctovl` | 0.0 | Minimum reciprocal overlap (0-1) |

WES-specific parameters can be set with `truvari_wes_*` prefix.

### Resource Limits

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_cpus` | 16 | Maximum CPUs per process |
| `max_memory` | 128.GB | Maximum memory per process |
| `max_time` | 240.h | Maximum time per process |

## Output Structure

```
results/
├── {run_name}/
│   ├── calls/                    # SV calls for each technology
│   │   ├── Illumina_WES/
│   │   │   └── sv/manta/
│   │   ├── Illumina_WGS/
│   │   │   └── sv/manta/
│   │   ├── PacBio/
│   │   │   ├── sv/cutesv/
│   │   │   └── sv/pbsv/
│   │   └── ONT/
│   │       ├── sv/cutesv/
│   │       └── sv/sniffles/
│   ├── benchmarking/             # Truvari results per target
│   │   ├── high_confidence/
│   │   ├── gene_panel/
│   │   └── wes_with_utr/
│   └── pipeline_info/           # Execution reports
│       ├── execution_report.html
│       ├── execution_timeline.html
│       └── execution_trace.txt
```

## Profiles

Available execution profiles:

| Profile | Description |
|---------|-------------|
| `singularity` | Use Singularity/Apptainer containers (recommended) |
| `docker` | Use Docker containers |
| `conda` | Use Conda environments |
| `test` | Run minimal test with small dataset |

Combine profiles with comma: `-profile test,docker`

## Advanced Usage

### Command-Line Parameter Override

```bash
nextflow run main.nf \
  -params-file params.yaml \
  --max_cpus 48 \
  --truvari_refdist 1000 \
  --outdir custom_output
```

### Running Specific Technologies

Omit unwanted BAM files from params file or set to empty string:

```yaml
illumina_wes_bam: /path/to/wes.bam    # Will run
illumina_wgs_bam: ''                   # Will skip
```

### Skip Specific Tools

```bash
nextflow run main.nf \
  -params-file params.yaml \
  --skip_pbsv true \              # Skip PBSV calling
  --skip_benchmarking true        # Skip Truvari benchmarking
```

### Target Region Simulation

For testing or exploring coverage, enable random target simulation:

```yaml
simulate_targets: true
num_simulated_targets: 10
simulated_target_size: 5000000
```

## Testing

Run the test profile to verify installation:

```bash
# Test with Docker
nextflow run main.nf -profile test,docker

# Test with Singularity
nextflow run main.nf -profile test,singularity
```

Expected runtime: ~5-10 minutes

## Troubleshooting

### BAM Index Not Found

Ensure `.bai` index files exist:
```bash
samtools index your_file.bam
```

### VCF Index Missing

Ensure VCF is bgzipped and indexed:
```bash
bgzip truth_set.vcf
tabix -p vcf truth_set.vcf.gz
```

### Container Mount Issues

For Singularity, use absolute paths and ensure auto-mounting is enabled in `nextflow.config`:
```groovy
singularity.autoMounts = true
```

### Memory Errors

Increase memory limits in command line or config:
```bash
nextflow run main.nf --max_memory 256.GB
```

## Container Images

The pipeline uses pre-built containers from:
- Docker Hub (for Docker profile)
- Locally cached images (for Singularity profile)

Container paths can be customized in `nextflow.config`:
```groovy
params {
    container_manta = "docker://biocontainers/manta:latest"
    container_cutesv = "docker://biocontainers/cutesv:latest"
    // ... etc
}
```

## Citation

If you use this pipeline, please cite:

- **Nextflow**: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology. https://doi.org/10.1038/nbt.3820
- **Manta**: Chen, X., et al. (2016). Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications. Bioinformatics.
- **CuteSV**: Jiang, T., et al. (2020). Long-read-based human genomic structural variation detection with cuteSV. Genome Biology.
- **Pbsv**: Pacific Biosciences. https://github.com/PacificBiosciences/pbsv
- **Sniffles**: Sedlazeck, F.J., et al. (2018). Accurate detection of complex structural variations using single-molecule sequencing. Nature Methods.
- **Truvari**: English, A.C., et al. (2022). Truvari: refined structural variant comparison preserves allelic diversity. Genome Biology.

## License

This pipeline is provided as-is for research purposes.

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request

## Support

For issues and questions:
- GitHub Issues: https://github.com/HudoGriz/SV_coding_regions_benchmark_nextflow/issues
