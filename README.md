# Benchmarking Structural Variants on Genome Intervals

A Nextflow DSL2 pipeline for systematic benchmarking of structural variant (SV) detection across multiple sequencing technologies and genomic interval categories using the Genome in a Bottle (GIAB) HG002 truth set.

This pipeline accompanies the manuscript:
> **How Diagnostic Target Selection Alters Structural Variant Benchmarking**
>
> Preprint: [10.21203/rs.3.rs-9179453/v1](https://doi.org/10.21203/rs.3.rs-9179453/v1)

## Overview

The pipeline evaluates SV detection performance across four sequencing platforms (Illumina WES, Illumina WGS, PacBio HiFi, ONT) within three genomic interval sets (high-confidence intervals, gene panel, exons+UTRs). Its simulation framework generates length- and chromosome-frequency-matched noncoding targets to measure the operational behavior of fragmented intervals. Because those targets do not match truth-SV composition or sequence context, they are an empirical interval null rather than a causal isolation of coding context.

### Pipeline Workflow

```
PREPARE_REFERENCES ─> SV_CALLING ─> BENCHMARKING ─┬─> SIMULATE_AND_BENCHMARK (optional)
                                                   └─> ANALYSIS_AND_PLOTS (optional)
```

1. **Prepare References** -- Validate and index reference genome, truth set, and target BED files
2. **SV Calling** -- Call structural variants with technology-appropriate callers
3. **Benchmarking** -- Compare calls against truth set using Truvari across all target intervals
4. **Simulation** -- Generate 500 random exon-like interval sets and benchmark against them
5. **Analysis** -- Compute statistics, percentile rankings, KDE outlier analysis, record-level target-transition evidence, and publication plots

### Supported Technologies and Callers

| Technology | SV Callers | Notes |
|-----------|-----------|-------|
| Illumina WES | Manta | Uses `--exome` flag; requires capture target BED |
| Illumina WGS | Manta | |
| PacBio HiFi | CuteSV, Pbsv | Pbsv can be skipped with `--skip_pbsv` |
| ONT | CuteSV, Sniffles | Sniffles supports tandem repeat annotation |

### Benchmarking Intervals

| Interval | Description |
|----------|-------------|
| **HCI** (High-Confidence Intervals) | GIAB-defined regions covering ~86% of the genome |
| **GP** (Gene Panel) | 3886 pediatric disorder genes |
| **EX+UTR** (Exons + UTRs) | GENCODE exonic regions with untranslated regions |

## Requirements

- Nextflow >= 23.04.0
- Container engine: Singularity/Apptainer (recommended) or Docker

## Quick Start

### 1. Prepare Data

The `preparation/` directory contains scripts to download all required GIAB data:

```bash
# Download data for both genome builds (~500 GB per build)
bash preparation/prepare.sh --genome all --outdir /path/to/data

# Or for a single build
bash preparation/prepare.sh --genome GRCh37 --outdir /path/to/data
```

This will:
- Download GIAB HG002 BAM files (Illumina WES/WGS, PacBio HiFi, ONT)
- Download reference genomes, truth sets, and annotations
- Create target BED files from GENCODE annotations
- Generate a ready-to-use `params_GRCh37.yaml` / `params_GRCh38.yaml`

### 2. Run the Pipeline

```bash
# Run with generated params file
nextflow run main.nf -params-file /path/to/data/GRCh37/params_GRCh37.yaml -profile singularity

# Resume if interrupted
nextflow run main.nf -params-file /path/to/data/GRCh37/params_GRCh37.yaml -profile singularity -resume
```

### 3. Custom Parameters File

If not using the preparation scripts, create a parameters file manually:

```yaml
# Required
fasta: /path/to/reference.fasta
benchmark_vcf: /path/to/truth_set.vcf.gz

# BAM files (provide only those you want to analyze)
illumina_wes_bam: /path/to/illumina_wes.bam
illumina_wgs_bam: /path/to/illumina_wgs.bam
pacbio_bam: /path/to/pacbio.bam
ont_bam: /path/to/ont.bam

# Target regions
high_confidence_targets: /path/to/high_conf.bed
gene_panel_targets: /path/to/gene_panel.bed
wes_utr_targets: /path/to/wes_utr.bed

# Optional
tandem_repeats: /path/to/tandem_repeats.bed
wes_sequencing_targets: /path/to/agilent_sureselect.bed.gz

# Output
outdir: ./results
run_name: my_benchmark
```

## Parameters

### Input Files

| Parameter | Required | Description |
|-----------|----------|-------------|
| `fasta` | Yes | Reference genome FASTA |
| `benchmark_vcf` | Yes | GIAB truth set VCF (.vcf.gz) |
| `illumina_wes_bam` | No | Illumina WES BAM file |
| `illumina_wgs_bam` | No | Illumina WGS BAM file |
| `pacbio_bam` | No | PacBio HiFi BAM file |
| `ont_bam` | No | Oxford Nanopore BAM file |
| `high_confidence_targets` | No | High-confidence regions BED |
| `gene_panel_targets` | No | Gene panel regions BED |
| `wes_utr_targets` | No | Exons + UTRs BED |
| `tandem_repeats` | No | Tandem repeat BED (improves Sniffles accuracy) |
| `wes_sequencing_targets` | No | WES capture targets BED.gz + .tbi (for Manta `--callRegions`) |

At least one BAM file must be provided.

### Pipeline Control

| Parameter | Default | Description |
|-----------|---------|-------------|
| `skip_benchmarking` | `false` | Skip Truvari benchmarking |
| `skip_pbsv` | `false` | Skip Pbsv caller for PacBio data |
| `simulate_targets` | `false` | Enable simulated interval analysis |
| `num_simulations` | `100` | Number of simulated interval sets to generate |
| `gather_statistics` | `false` | Generate publication plots and statistics tables |
| `generate_transition_evidence` | `false` | Trace HCI-to-target outcome changes through exact Truvari MatchIds |

### Truvari Parameters

Default parameters for SV comparison. Separate `truvari_wes_*` parameters allow different thresholds for WES data.

| Parameter | Default | WES Default | Description |
|-----------|---------|-------------|-------------|
| `truvari_refdist` | 500 | 500 | Max reference distance (bp) |
| `truvari_pctsize` | 0.7 | 0.7 | Min size similarity (0-1) |
| `truvari_pctseq` | 0.0 | 0.0 | Min sequence similarity (0-1) |
| `truvari_pctovl` | 0.0 | 0.0 | Min reciprocal overlap (0-1) |

All Truvari runs include `--bench-overlaps 1 --passonly --dup-to-ins` flags. The pipeline uses a [modified Truvari](https://github.com/CISLD/truvari) that allows partial overlap with target intervals. The value of `--bench-overlaps` is the minimum number of positions a call must share with a target interval; `1` is the one-base intersection the published results use, and `0` restores stock containment.

### Resource Limits

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_cpus` | 24 | Maximum CPUs per process |
| `max_memory` | 128.GB | Maximum memory per process |
| `max_time` | 48.h | Maximum time per process |

## Output Structure

```
{outdir}/
├── sv_calls/                        # SV caller output VCFs
│   ├── Illumina_WES/Manta/
│   ├── Illumina_WGS/Manta/
│   ├── PacBio/
│   │   ├── CuteSV/
│   │   └── PBSV/
│   └── ONT/
│       ├── CuteSV/
│       └── Sniffles/
├── real_intervals/                  # Truvari benchmarks on real target sets
│   └── {technology}-{caller}-{target}/
├── simulations/                     # Simulated interval analysis (if enabled)
│   ├── simulated_targets/           # Generated BED files
│   └── benchmarks/                  # Truvari results per simulation
├── statistics/                      # Plots and tables (if enabled)
│   ├── plots/
│   │   ├── bar_plot.png
│   │   ├── bar_plot_sim_diff.png
│   │   └── facets_plot.png
│   └── tables/
│       ├── truvari_metrics_real_intervals.tsv
│       ├── truvari_metrics_simulated_intervals.tsv
│       └── truvari_metrics_simulated_intervals_raw.tsv
└── pipeline_info/                   # Nextflow execution reports
    ├── execution_report.html
    ├── execution_timeline.html
    └── execution_trace.txt
```

## Profiles

| Profile | Description |
|---------|-------------|
| `singularity` | Singularity/Apptainer containers (recommended for HPC) |
| `docker` | Docker containers |
| `test` | Local test with small dataset (requires `test_data/` directory) |
| `test_nfcore` | Remote nf-core test data for CI (SV calling only, no benchmarking) |

Combine profiles: `-profile singularity` or `-profile test_nfcore,docker`

### Target-transition evidence and publication figures

The analysis image is published, so nothing needs building; Nextflow pulls
`library://blazv/benchmark-sv/python-r-analysis:py3.11-r4.4.1` on first use. It
contains the record-level MatchId audit, batch merger, factor attribution,
occurrence mapping, and publication plotting dependencies. Enable the integrated
workflow with:

```bash
nextflow run . -profile singularity \
    --generate_transition_evidence \
    --reference_assembly GRCh37
```

To iterate on the image, rebuild it from the repository root and point the
pipeline at the local file:

```bash
bin/build_python_r_analysis_container.sh
nextflow run . -profile singularity \
    --generate_transition_evidence \
    --reference_assembly GRCh37 \
    --analysis_container containers/python-r-analysis.sif
```

Evidence is published under `target_transition_evidence/` in the run output.
The default comparison is `high_confidence` versus `wes_utr`, labelled
`EX+UTR` in publication tables. Override these with
`--transition_hci_target`, `--transition_target`, and
`--transition_target_label`.

Nextflow itself runs on the host. On the KISLD HPC installation used for the
paper, load it with `module load anaconda` followed by `conda activate nf-core`.
Callers, Truvari, and the combined Python/R analysis environment remain separate
process containers; the workflow is not launched from inside a container.

> **Container engine.** The two custom images are published to a Singularity
> library, and `library://` is not an OCI registry. Benchmarking and analysis
> therefore require Singularity or Apptainer. The `docker`, `podman` and
> `charliecloud` profiles can run SV calling, whose images all come from
> `quay.io`, but not the Truvari or analysis stages.

### Running the drivers elsewhere

The scripts in `bin/` contain no absolute paths. Everything site-specific is an
environment variable, so they run unmodified on any machine:

| Variable | Default | Purpose |
|----------|---------|---------|
| `SV_DATA_ROOT` | *required* | Directory holding the prepared per-assembly data |
| `SV_PROFILE` | `singularity` | Nextflow profile to run with |
| `SV_HPC_CONFIG` | unset | Extra Nextflow config for a local cluster |
| `SV_ENV_MODULE` | unset | Environment module to load before running |
| `SV_CONDA_ENV` | unset | Conda environment to activate before running |
| `TRUVARI_SIF` | published image | Local `.sif` or registry URI overriding the Truvari container |
| `ANALYSIS_SIF` | published image | Local `.sif` overriding the analysis container |
| `SV_IMAGE_CACHE` | inside the output directory | Where pulled images are cached |

Leaving `TRUVARI_SIF` and `ANALYSIS_SIF` unset is the reproducible choice: the
pipeline then uses the published, immutable tags in `nextflow.config`. Setting
either one is recorded in the run manifest along with its checksum.

To reproduce the published runs on the machine that produced them:

```bash
export SV_DATA_ROOT=/path/to/prepared/data
export SV_ENV_MODULE=anaconda SV_CONDA_ENV=nf-core SV_PROFILE=cpu
export SV_HPC_CONFIG=/home/Software/configs/nextflow_local/configs/conf/kisld_hpc.config
bin/run_clean_dated_benchmark.sh 2026-08-03
```

## Repository Structure

```
main.nf                         # Pipeline entry point
nextflow.config                 # Main configuration
nextflow_schema.json            # JSON Schema for parameter validation
conf/
  modules.config                # Per-process containers, publishDir, ext.args
  test.config                   # Local test profile
  test_nfcore.config            # Remote CI test profile
workflows/
  prepare_references.nf         # Reference/index validation
  sv_calling.nf                 # SV caller orchestration
  benchmarking.nf               # Truvari benchmarking across intervals
  simulate_and_benchmark.nf     # Simulated interval generation + benchmarking
  target_transition_evidence.nf # Record-level HCI-to-target outcome audit
  analysis_and_plots.nf         # Statistics and plot generation
modules/
  local/
    simulate_targets.nf         # Random exon-like interval simulation
    gather_statistics.nf        # R-based statistics and plotting
    truvari_refine.nf           # Truvari refine on benchmark output
    target_transition_evidence.nf # Audit, merge, and plot processes
  nf-core/                      # Pinned nf-core modules
containers/
  Singularity.python-r-analysis # Combined Python/R analysis image definition
bin/R/
  simulate_targets.R            # Simulation algorithm (GenomicRanges-based)
  paper_plots.R                 # Publication plot generation
  functions.R                   # Shared R utilities
bin/python/
  target_transition_audit.py    # MatchId-based HCI-to-target transition audit
  merge_transition_audits.py    # Merge batched audits into evidence tables
  plot_target_boundary_mechanisms.py # Boundary-mechanism figure
  exutr_factor_attribution.py   # EX+UTR factor attribution
  sv_occurrence_map.py          # Per-record occurrence mapping
  pad_target_bed.py             # Pad, clip, and merge a target BED
  padding_transition_recovery.py # Trace losses across padded targets
  summarize_padding_sensitivity.py # Padding sensitivity tables and plots
  compare_metric_tables.py      # Diff two gather-statistics tables
  target_composition.py         # Type and size composition of scored truth
  update_submission_tables.py   # Refresh supplementary tables from a run
preparation/
  prepare.sh                    # Master data download wrapper
  download_and_prep_GRCh37.sh   # GRCh37 data acquisition
  download_and_prep_GRCh38.sh   # GRCh38 data acquisition
  generate_params.sh            # Auto-generate params YAML from downloaded data
  create_gencode_target_bed.R   # Create exon+UTR BED from GENCODE GTF
```

## Testing

```bash
# CI test with remote nf-core data (Docker)
nextflow run main.nf -profile test_nfcore,docker --outdir test_results

# Local test (requires test_data/ directory)
nextflow run main.nf -profile test,singularity --outdir test_results
```

## Citation

If you use this pipeline, please cite:

- **Truvari**: English, A.C., et al. (2022). Truvari: refined structural variant comparison preserves allelic diversity. *Genome Biology*, 23, 271.
- **Nextflow**: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*, 35, 316-319.
- **Manta**: Chen, X., et al. (2016). Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications. *Bioinformatics*, 32, 1220-1222.
- **CuteSV**: Jiang, T., et al. (2020). Long-read-based human genomic structural variation detection with cuteSV. *Genome Biology*, 21, 189.
- **Pbsv**: Pacific Biosciences. https://github.com/PacificBiosciences/pbsv
- **Sniffles**: Smolka, M., et al. (2024). Detection of mosaic and population-level structural variants with Sniffles2. *Nature Biotechnology*, 42, 1571-1580.
- **GIAB**: Zook, J.M., et al. (2020). A robust benchmark for detection of germline large deletions and insertions. *Nature Biotechnology*, 38, 1347-1355.

## License

This pipeline is provided as-is for research purposes.
