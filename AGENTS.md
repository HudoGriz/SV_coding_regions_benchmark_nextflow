# AGENTS.md — SV Calling & Benchmarking Nextflow Pipeline

## Project Overview

Nextflow DSL2 pipeline for structural variant (SV) calling and benchmarking across
multiple sequencing technologies (Illumina WES/WGS, PacBio HiFi, ONT), built
following nf-core conventions. Uses Genome in a Bottle (GIAB) truth sets for
benchmarking with Truvari. Requires Nextflow >= 23.04.0.

## Repository Layout

```
main.nf                    # Pipeline entry point
nextflow.config            # Main config (params, resources, profiles)
nextflow_schema.json       # JSON Schema for pipeline parameters
conf/
  modules.config           # Per-process config (containers, publishDir, ext.args)
  test.config              # Minimal local test profile
  test_nfcore.config       # Remote nf-core test data profile
workflows/                 # Sub-workflows (SCREAMING_SNAKE_CASE names)
modules/
  local/                   # Custom process definitions (simulate_targets, gather_statistics,
                           #   truvari_refine, target_transition_evidence)
  nf-core/                 # Pinned nf-core modules (do NOT edit directly)
lib/WorkflowHelp.groovy    # Groovy helper class
bin/R/                     # R scripts (simulate_targets.R, paper_plots.R, functions.R)
bin/python/                # Analysis scripts shipped in the analysis container
bin/*.sh                   # Container build and evidence-regeneration drivers
containers/
  Singularity.python-r-analysis  # Combined Python/R analysis image definition
preparation/               # Shell scripts for data download/prep
```

## Build / Run / Test Commands

```bash
# On this HPC host, enable Nextflow first
module load anaconda
conda activate nf-core

# Quick validation / syntax check (the CI lint step, no containers needed)
nextflow run . -profile test_nfcore --help

# Full integration test with nf-core remote test data (Docker)
nextflow run . -profile test_nfcore,docker --outdir test_results -resume

# Full integration test with local test data (requires test_data/ directory)
nextflow run . -profile test,docker --outdir test_results

# Run a single nf-core module test (requires nf-test installed)
nf-test test modules/nf-core/truvari/bench/tests/main.nf.test

# View resolved configuration for a profile
nextflow config -profile test_nfcore,docker

# Clean Nextflow work directory
rm -rf work/ .nextflow/ .nextflow.log*
```

**CI** (`.github/workflows/ci.yml`) runs 5 jobs on push/PR to main/dev/master:
lint (syntax check), nf-core-lint (optional, non-blocking), module-check,
profile-check, and run-test (full pipeline with `test_nfcore,docker`).

No unit tests exist for the pipeline itself — testing is integration-based.
Module tests use nf-test: `modules/nf-core/<tool>/tests/main.nf.test`.

## Linting & Formatting

No formal linter or formatter is enforced. The CI optionally runs `nf-core pipelines lint`
(continue-on-error). Follow the conventions observed in existing code:

- **Indentation**: 4 spaces (Nextflow, Groovy, R, shell). No tabs.
- **Line length**: No hard limit; keep lines readable (~120 chars).
- **Trailing whitespace**: Avoid. **File encoding**: UTF-8, LF line endings.

## Code Style — Nextflow / Groovy

### Naming Conventions

| Element | Convention | Example |
|---------|-----------|---------|
| Workflows / Processes | `SCREAMING_SNAKE_CASE` | `SV_CALLING`, `TRUVARI_BENCH` |
| Channels | `ch_` prefix + `snake_case` | `ch_all_vcfs`, `ch_fasta_fai` |
| Parameters | `snake_case` | `params.illumina_wes_bam`, `params.skip_pbsv` |
| Local variables | `snake_case` | `def is_remote`, `def pacbio_bam` |
| Meta map keys | `snake_case` | `meta.technology`, `meta.tool`, `meta.target` |
| Module aliases | `SCREAMING_SNAKE_CASE` | `MANTA_GERMLINE as MANTA_WES` |
| Groovy classes | `PascalCase` | `WorkflowHelp` |

### Imports

- `include { ... } from '...'` at the top of each file, after the header comment.
- Group: nf-core modules first, then local modules.
- Use `as` aliases when the same module is used multiple times:
  ```groovy
  include { CUTESV as CUTESV_PACBIO } from '../modules/nf-core/cutesv/main'
  include { CUTESV as CUTESV_ONT }    from '../modules/nf-core/cutesv/main'
  ```

### Channel Patterns

- Empty channels: `Channel.empty()`. Single-value: `Channel.value()`. Lists: `Channel.from()`.
- File channels: `Channel.fromPath(path, checkIfExists: true)`.
- Combine with `.combine()`, merge with `.mix()`, transform with `.map { }`.

### Process Definitions (local modules)

```groovy
process PROCESS_NAME {
    tag "descriptive tag"
    label 'process_low'        // process_low | process_medium | process_high

    input:
    val some_value
    path some_file

    output:
    path "output/*", emit: named_output

    script:
    """
    bash_commands_here
    """
}
```

- Process config (containers, publishDir, ext.args) goes in `conf/modules.config`,
  NOT in the process definition. Use `withName:` selectors.
- Keep process scripts minimal; complex logic goes in `bin/` scripts.

### Section Comments

Use banner style for major sections:
```groovy
/*
========================================================================================
    SECTION TITLE
========================================================================================
*/
```
Use `//` for inline comments. Document channel shapes on `take:`/`emit:` lines
(e.g., `// channel: [meta, vcf, tbi]`).

### Error Handling

- Fatal errors: `error("message")` or `log.error` + `error()` for detailed messages.
- Non-fatal warnings: `log.warn`.
- Validate required parameters early in `main.nf` before calling sub-workflows.
- Retry strategy in `nextflow.config`: `errorStrategy = { task.exitStatus in [143,137,104,134,139] ? 'retry' : 'finish' }`.

### Meta Maps

Every VCF channel carries a `meta` map as the first tuple element:
```groovy
[id: 'sample_id', technology: 'PacBio', tool: 'CuteSV', target: 'high_confidence']
```
Always propagate `meta` through transformations. Add keys via `meta + [key: value]`.

## Code Style — R Scripts (`bin/R/`)

- `library()` imports at top. Packages: `ggplot2`, `jsonlite`, `GenomicRanges`, `rtracklayer`, `parallel`.
- `commandArgs(trailingOnly = TRUE)` for CLI args. Use `<-` for assignment. `snake_case` everywhere.

## nf-core Modules

Pinned in `modules.json` at specific git SHAs. **Do not edit directly.** Update with:
```bash
nf-core modules update <module_name>
```
Override behavior in `conf/modules.config` via `withName:` blocks.

## Container Configuration

- Default registry: `quay.io` (Docker, Singularity, Apptainer, Podman, Charliecloud).
- Custom Singularity containers from `library://blazv/benchmark-sv/`:
  `r-env:4-4-1` (R environment), `truvari_modded:latest` (modified Truvari).
- Analysis processes resolve their image from `params.analysis_container`. Build it with
  `bin/build_python_r_analysis_container.sh`; the definition is
  `containers/Singularity.python-r-analysis`. Built `.sif` files are gitignored.
- `conf/modules.config` is included after the `params` block in `nextflow.config` so that
  `container = { params.analysis_container }` closures resolve without undefined-parameter warnings.
- Per-process container overrides in `conf/modules.config`.

## Key Parameters

Defined in `nextflow.config` with defaults, documented in `nextflow_schema.json`.
Boolean flags: `skip_benchmarking`, `skip_pbsv`, `simulate_targets`, `gather_statistics`.
Use `null` as default for optional file paths.
