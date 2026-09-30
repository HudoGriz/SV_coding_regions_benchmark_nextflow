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
                                                   ├─> SENSITIVITY_BENCHMARKS (optional)
                                                   ├─> TARGET_TRANSITION_EVIDENCE (optional)
                                                   └─> ANALYSIS_AND_PLOTS (optional)
```

1. **Prepare References** -- Validate and index reference genome, truth set, and target BED files
2. **SV Calling** -- Call structural variants with technology-appropriate callers
3. **Benchmarking** -- Compare calls against truth set using Truvari across all target intervals
4. **Simulation** -- Generate 500 random exon-like interval sets and benchmark against them
5. **Sensitivity** -- Re-score the real targets with one setting changed at a time: matching thresholds, full containment, candidate-only `--extend`, and symmetric padding
6. **Analysis** -- Compute statistics, percentile rankings, KDE outlier analysis, record-level target-transition evidence, and publication plots

### Supported Technologies and Callers

| Technology | SV Callers | Notes |
|-----------|-----------|-------|
| Illumina WES | Manta | Uses `--exome` flag; requires capture target BED |
| Illumina WGS | Manta, Delly | Delly (v1.7.3) runs with its exclude template; skip with `--skip_delly` |
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
| `skip_delly` | `false` | Skip Delly on Illumina WGS |
| `delly_exclude` | `null` | Delly exclude template (downloaded by the preparation scripts) |
| `sensitivity_benchmarks` | `false` | Re-score the real targets under alternative settings (`sensitivity_*` parameters) |
| `simulate_targets` | `false` | Enable simulated interval analysis |
| `num_simulations` | `100` | Number of simulated interval sets to generate |
| `gather_statistics` | `false` | Generate publication plots and statistics tables |
| `generate_transition_evidence` | `false` | Trace HCI-to-target outcome changes through exact Truvari MatchIds |

### Truvari Parameters

Default parameters for SV comparison. Separate `truvari_wes_*` parameters allow different thresholds for WES data.

| Parameter | Default | WES Default | Description |
|-----------|---------|-------------|-------------|
| `truvari_refdist` | 500 | 500 | Max reference distance (bp) |
| `truvari_pctsize` | 0.7 | 0.7 (generated params files: 0) | Min size similarity (0-1) |
| `truvari_pctseq` | 0.0 | 0.0 | Min sequence similarity (0-1) |
| `truvari_pctovl` | 0.0 | 0.0 | Min reciprocal overlap (0-1) |

All Truvari runs include `--bench-overlaps 1 --bnddist -1 --passonly --dup-to-ins`, and keep
Truvari's size defaults: truth records must be at least `--sizemin` 50 bp, candidates at least
`--sizefilt` 30 bp, and both at most `--sizemax` 50 kb. A candidate of 30-49 bp may match a truth
record but is dropped, not counted as a false positive, when it does not. `--refdist` is satisfied
when the candidate's span lies within that distance of the truth span; `--pctseq` is applied only
when both records are sequence-resolved. Inversions are not converted and are scored as their own
type. Every benchmark directory keeps its full configuration in `<prefix>/params.json`.

Before any benchmark, calls genotyped homozygous reference (`0/0`, `0|0`, haploid `0`) are removed (`exclude_homref_calls`, default `true`): the caller is stating that the sample does not carry them, and the truth sets count only records that carry an ALT allele. Calls without a genotype (`./.`, all cuteSV calls) are kept, which is why Truvari's own `--no-ref` is not used. A caller with no such call is benchmarked on its original VCF. The scored VCFs and the per-caller counts are in `benchmarked_calls/`. `--exclude_homref_calls false` scores them like any other call. The pipeline uses a [modified Truvari](https://github.com/CISLD/truvari) that allows partial overlap with target intervals. The value of `--bench-overlaps` is the minimum number of positions a call must share with a target interval; `1` is the one-base intersection the published results use, and `0` restores stock containment.

### Resource Limits

| Parameter | Default | Description |
|-----------|---------|-------------|
| `max_cpus` | 24 | Maximum CPUs per process |
| `max_memory` | 128.GB | Maximum memory per process |
| `max_time` | 48.h | Maximum time per process |

## Output Structure

```
{outdir}/
├── benchmarked_calls/               # VCFs as scored (0/0 calls removed) and the removal counts
├── sv_calls/                        # SV caller output VCFs
│   ├── Illumina_WES/Manta/
│   ├── Illumina_WGS/Manta/
│   ├── Illumina_WGS/Delly/
│   ├── PacBio/
│   │   ├── CuteSV/
│   │   └── PBSV/
│   └── ONT/
│       ├── CuteSV/
│       └── Sniffles/
├── real_intervals/                  # Truvari benchmarks on real target sets
│   └── {technology}-{caller}-{target}/
├── sensitivity/                     # Sensitivity benchmarks (if enabled)
│   ├── {setting}/{technology}/{caller}/{target}/   # refdist*, pctsize*, pctseq*, containment, extend*, pad*
│   └── target_beds/                 # Padded boundary-target BEDs
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

Nextflow (25.04 or later) runs on the host and must be on `PATH`; how it gets
there is up to the installation. Callers, Truvari, and the combined Python/R analysis environment remain separate
process containers; the workflow is not launched from inside a container.

> **Container engine.** The two custom images are published to a Singularity
> library, and `library://` is not an OCI registry. Benchmarking and analysis
> therefore require Singularity or Apptainer. The `docker`, `podman` and
> `charliecloud` profiles can run SV calling, whose images all come from
> `quay.io`, but not the Truvari or analysis stages.

### Running the study

One command runs every analysis of the study for one assembly: SV calling
(including Delly), benchmarking on HCI, GP and EX+UTR, the 500 simulated interval
sets, statistics and plots, target-transition evidence, the sensitivity
benchmarks and the post-hoc analyses. The `study` profile sets those options;
each one remains an ordinary parameter that can be switched off.

```bash
bash preparation/generate_params.sh --genome GRCh38 --datadir /path/to/GRCh38/data \
    --outfile params_GRCh38.yaml
nextflow run . -params-file params_GRCh38.yaml -profile singularity,study \
    --reference_assembly GRCh38 --outdir results-GRCh38
```

The post-hoc analyses (`--posthoc_analyses`, on in `study`) run as pipeline
tasks after the benchmarks and publish under `<outdir>/posthoc/`:

- every benchmark's metrics and Truvari parameters in one table, optionally
  compared row by row with another run (`--posthoc_compare_results <dir>`);
- SV-type accounting from the scored VCF to the HCI benchmark;
- the recall and precision decomposition (composition versus label transitions)
  with the candidate-side audit of false positives, per pipeline;
- post-matching stratification, checked against the composition-only values;
- block-bootstrap intervals and rarefied percentile ranks;
- transition audits across the threshold grid and recovery under `--extend` and
  padding;
- composition standardisation of the simulated metrics;
- simulation fidelity (GIAB v3.3 stratifications by default;
  `--posthoc_segdups` and `--posthoc_lowmappability` override them).

They need a truth set and `--simulate_targets true`; with the sensitivity
benchmarks they also need `--generate_transition_evidence true`.

The SVanalyzer second-comparator check is deliberately not part of the pipeline.
Run it on a finished run with `bin/run_svanalyzer_posthoc.sh <run_root> <assembly>`
(`SVANALYZER_SIF` built from `containers/Singularity.svanalyzer`).

GRCh38 BAMs are prepared with `preparation/build_grch38_analysis_bams.sh`, which
restricts them to the contigs of the analysis reference and filters nothing
else: discordant pairs, supplementary alignments and reads with an unmapped mate
all stay. A second stage (`strip_absent_sa_entries.py`) removes the SA-tag
entries that still point at a dropped contig, which pbsv otherwise aborts on; it
keeps every record. A BAM whose header already matches the reference (ONT) is
used as distributed. The `*.analysis_contigs.sa_filtered.bam` files in
`data/analysis_bams/` are the pipeline inputs.

### Run driver

`bin/run_benchmark.sh <label> <assembly> [nextflow arguments ...]` wraps the
command above for dated, provenance-recorded runs. Each assembly runs from its
own directory under `$SV_DATA_ROOT/<label>/`; the driver refuses to write into
existing results unless `SV_RESUME=1`, and records the params file, the
execution config, the image checksums and a hash of the pipeline files, so a
run can be matched to a commit even where the execution host has no git. Extra
arguments go to Nextflow unchanged.

| Variable | Default | Purpose |
|----------|---------|---------|
| `SV_DATA_ROOT` | *required* | Directory holding the prepared per-assembly data |
| `SV_PARAMS_FILE` | `$SV_DATA_ROOT/<assembly>/params_<assembly>.yaml` | Params file |
| `SV_PROFILE` | `singularity,study` | Nextflow profile(s) |
| `SV_HPC_CONFIG` | unset | Extra Nextflow config for the execution environment (executor, queue, container cache) |
| `SV_ENV_MODULE` | unset | Environment module to load before running |
| `SV_CONDA_ENV` | unset | Conda environment to activate before running |
| `TRUVARI_SIF` | published image | Local `.sif` or registry URI overriding the Truvari container |
| `ANALYSIS_SIF` | published image | Local `.sif` overriding the analysis container |
| `SV_RESUME` | unset | `1` resumes the run instead of refusing |

Leaving `TRUVARI_SIF` and `ANALYSIS_SIF` unset is the reproducible choice: the
pipeline then uses the published, immutable tags in `nextflow.config`. Setting
either one is recorded in the run manifest along with its checksum.

```bash
export SV_DATA_ROOT=/path/to/prepared/data
export SV_HPC_CONFIG=/path/to/your/site.config   # optional
for asm in GRCh37 GRCh38; do
    bin/run_benchmark.sh "$(date +%F)" "$asm"
done
```

The driver assumes no scheduler. To run it as a batch job, submit it with your
scheduler's own command; launchers specific to one site belong in the
git-ignored `local/` directory, not in the repository.

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
- **Delly**: Rausch, T., et al. (2012). DELLY: structural variant discovery by integrated paired-end and split-read analysis. *Bioinformatics*, 28, i333-i339.
- **Manta**: Chen, X., et al. (2016). Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications. *Bioinformatics*, 32, 1220-1222.
- **CuteSV**: Jiang, T., et al. (2020). Long-read-based human genomic structural variation detection with cuteSV. *Genome Biology*, 21, 189.
- **Pbsv**: Pacific Biosciences. https://github.com/PacificBiosciences/pbsv
- **Sniffles**: Smolka, M., et al. (2024). Detection of mosaic and population-level structural variants with Sniffles2. *Nature Biotechnology*, 42, 1571-1580.
- **GIAB**: Zook, J.M., et al. (2020). A robust benchmark for detection of germline large deletions and insertions. *Nature Biotechnology*, 38, 1347-1355.

## License

This pipeline is provided as-is for research purposes.
