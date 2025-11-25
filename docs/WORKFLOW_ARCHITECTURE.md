# Workflow Architecture

## Overview

The pipeline is organized into **6 independent sub-workflows**, each with a specific purpose. This modular design makes the code easier to maintain, test, and extend.

## Workflow Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                         MAIN WORKFLOW                            │
│                           (main.nf)                              │
└──────────────────────────┬──────────────────────────────────────┘
                           │
        ┌──────────────────┼──────────────────┬──────────────────┐
        │                  │                  │                  │
        ▼                  ▼                  ▼                  ▼
┌───────────────┐  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐
│  PREPARE_     │  │  SV_CALLING  │  │ BENCHMARKING │  │  SIMULATE_   │
│  REFERENCES   │  │              │  │              │  │  AND_        │
│               │  │              │  │              │  │  BENCHMARK   │
└───────┬───────┘  └──────┬───────┘  └──────┬───────┘  └──────┬───────┘
        │                 │                 │                 │
        │                 │                 │                 │
        │                 │                 │                 │
┌───────▼─────────────────▼─────────────────▼─────────────────▼───────┐
│                                                                       │
│                      ANALYSIS_AND_PLOTS                               │
│                   (Optional aggregation)                              │
│                                                                       │
└───────────────────────────────────────────────────────────────────────┘
```

## Sub-Workflow Details

### 1. PREPARE_REFERENCES (`workflows/prepare_references.nf`)

**Purpose**: Load and prepare all reference files

**Inputs** (from params):
- `params.fasta` - Reference genome FASTA
- `params.benchmark_vcf` - Truth VCF (optional)
- `params.high_confidence_targets` - BED file
- `params.gene_panel_targets` - BED file
- `params.wes_utr_targets` - BED file
- `params.tandem_repeats` - BED file (optional)

**Outputs**:
- `fasta` - Reference FASTA channel
- `fasta_fai` - FAI index channel (auto-created if missing)
- `benchmark_vcf` - Truth VCF channel
- `benchmark_vcf_tbi` - Truth VCF index channel
- `targets` - Channel: [target_name, bed_file]
- `tandem_repeats` - Tandem repeats channel

**Key Features**:
- Automatic FAI index creation using SAMTOOLS_FAIDX
- Remote file support (HTTP/HTTPS/FTP)
- Conditional benchmark VCF handling
- Clean channel organization for downstream workflows

---

### 2. SV_CALLING (`workflows/sv_calling.nf`)

**Purpose**: Call structural variants across multiple sequencing technologies

**Inputs**:
- `ch_fasta` - Reference FASTA
- `ch_fasta_fai` - Reference FAI index
- `ch_tandem_repeats` - Tandem repeats BED (for Sniffles)

**Technology Support**:

| Technology | Tools | Notes |
|------------|-------|-------|
| Illumina WES | Manta | WES-specific parameters |
| Illumina WGS | Manta | Standard parameters |
| PacBio | CuteSV + PBSV | PBSV optional (`--skip_pbsv`) |
| ONT | CuteSV + Sniffles | Sniffles uses tandem repeats |

**Outputs**:
- `vcfs` - Channel: [meta, vcf, tbi] for all SV calls

**Key Features**:
- Conditional execution based on input BAMs
- Automatic VCF compression and indexing
- Technology and tool metadata propagation
- Parallel execution across all technologies

---

### 3. BENCHMARKING (`workflows/benchmarking.nf`)

**Purpose**: Benchmark SV calls against truth set using Truvari

**Inputs**:
- `ch_vcfs` - SV caller VCF files (from SV_CALLING)
- `ch_benchmark_vcf` - Truth VCF
- `ch_benchmark_vcf_tbi` - Truth VCF index
- `ch_targets` - Target region sets
- `ch_fasta` - Reference FASTA
- `ch_fasta_fai` - Reference FAI index

**Benchmarking Strategy**:
- Runs Truvari on **predefined target sets**:
  - `high_confidence` - High-confidence regions
  - `gene_panel` - Gene panel regions
  - `wes_utr` - Whole exome + UTR regions

**Outputs**:
- `summary` - Truvari summary.json files for each combination of:
  - Technology (Illumina_WES, Illumina_WGS, PacBio, ONT)
  - Tool (Manta, CuteSV, PBSV, Sniffles)
  - Target set (high_confidence, gene_panel, wes_utr)

**Key Features**:
- Automatic parameter selection (WES vs WGS)
- Technology-specific Truvari parameters
- Per-target, per-tool benchmarking
- Metadata tracking through pipeline

---

### 4. SIMULATE_AND_BENCHMARK (`workflows/simulate_and_benchmark.nf`)

**Purpose**: Create random target regions and benchmark against them

**⚠️ This is SEPARATE from the main BENCHMARKING workflow**

**Inputs**:
- `ch_fasta` - Reference FASTA
- `ch_fasta_fai` - Reference FAI index
- `ch_gencode_gtf` - GENCODE GTF annotation
- `ch_benchmark_vcf` - Truth VCF
- `ch_benchmark_vcf_tbi` - Truth VCF index
- `ch_vcfs` - SV caller VCF files (from SV_CALLING)
- `num_simulations` - Number of simulations (default: 100)

**Workflow**:
1. **SIMULATE_TARGETS**: Extracts genes from GTF and creates random combinations
2. **TRUVARI_BENCH**: Benchmarks each SV caller against each simulated target

**Outputs**:
- `simulated_beds` - Channel of simulated BED files
- `truvari_results` - Truvari summary.json files for simulations

**Key Features**:
- Creates diverse target sets similar to exome+UTR
- Each simulation has different gene combinations
- Enables statistical analysis of performance variation
- Useful for understanding sensitivity to target selection

**Differences from main BENCHMARKING**:
| Feature | BENCHMARKING | SIMULATE_AND_BENCHMARK |
|---------|--------------|------------------------|
| Target sets | Predefined (3 fixed sets) | Random (100+ simulations) |
| Execution | Always (if benchmark_vcf provided) | Optional (`--simulate_targets`) |
| Purpose | Evaluate on known targets | Statistical robustness testing |
| Parameters | Uses standard/WES-specific | Uses standard parameters |

---

### 5. ANALYSIS_AND_PLOTS (`workflows/analysis_and_plots.nf`)

**Purpose**: Aggregate results and generate visualizations

**Inputs**:
- `ch_truvari_results` - All Truvari summary files
- `ch_run_dir` - Output directory

**Outputs**:
- Aggregated statistics CSV
- Performance plots
- Comparative analyses

**Key Features**:
- Combines results from main benchmarking and simulations
- Generates publication-ready plots
- Calculates aggregate metrics

---

### 6. PREPARE_GIAB_RESOURCES (`workflows/prepare_giab_resources.nf`)

**Purpose**: Download GIAB truth sets and annotations

**Mode**: Data preparation (optional)

**Outputs**:
- GIAB truth VCF files
- Annotation files (tandem repeats, etc.)
- Target BED files

---

## Execution Flow

### Standard Analysis Run

```
1. PREPARE_REFERENCES
   └─> Loads files, creates indices
       │
2. SV_CALLING
   └─> Calls SVs across all technologies in parallel
       │
3. BENCHMARKING (if --benchmark_vcf provided)
   └─> Benchmarks against predefined targets
       │
4. SIMULATE_AND_BENCHMARK (if --simulate_targets enabled)
   └─> Creates random targets and benchmarks
       │
5. ANALYSIS_AND_PLOTS (if --gather_statistics enabled)
   └─> Aggregates and visualizes results
```

### Data Preparation Run

```
1. PREPARE_GIAB_RESOURCES (if --prepare_giab_resources)
   └─> Downloads minimal GIAB data

OR

1. PREPARE_DATA_COMPLETE_GRCH37/38 (if --prepare_complete_data)
   └─> Downloads complete dataset including BAMs
```

## Parallelization

The pipeline automatically parallelizes at multiple levels:

1. **Technology-level**: All technologies run in parallel
   - Illumina WES
   - Illumina WGS
   - PacBio (CuteSV + PBSV in series)
   - ONT (CuteSV + Sniffles in series)

2. **Benchmarking-level**: All VCF × Target combinations run in parallel
   - Example: 6 VCFs × 3 targets = 18 parallel Truvari jobs

3. **Simulation-level**: All simulations run in parallel
   - Example: 100 simulations × 6 VCFs = 600 parallel jobs

## Adding New Components

### To add a new SV caller:

1. Add nf-core module to `modules/nf-core/`
2. Edit `workflows/sv_calling.nf`:
   ```groovy
   if (params.new_tool_bam) {
       NEW_TOOL(
           ch_input,
           ch_fasta,
           ch_fasta_fai
       )
       ch_all_vcfs = ch_all_vcfs.mix(NEW_TOOL.out.vcf)
   }
   ```

### To add a new benchmarking strategy:

1. Create new workflow: `workflows/my_benchmark.nf`
2. Include in `main.nf`:
   ```groovy
   include { MY_BENCHMARK } from './workflows/my_benchmark'
   
   MY_BENCHMARK(
       SV_CALLING.out.vcfs,
       ch_fasta,
       ch_fasta_fai
   )
   ```

### To modify target selection:

- **Fixed targets**: Edit `workflows/benchmarking.nf`
- **Simulated targets**: Edit `bin/simulate_targets.py`

## Summary

The modular architecture provides:

✅ **Clear separation of concerns**: Each workflow has one job  
✅ **Independent testing**: Test workflows in isolation  
✅ **Easy extension**: Add features without touching main workflow  
✅ **Parallel execution**: Automatic parallelization at multiple levels  
✅ **Maintainability**: Changes isolated to specific components  
✅ **Reusability**: Sub-workflows can be imported by other pipelines  

The simulation workflow (`SIMULATE_AND_BENCHMARK`) is **already separate** from the main benchmarking workflow (`BENCHMARKING`), giving you flexibility to:
- Run only on predefined targets (fast, focused)
- Run only simulations (statistical robustness)
- Run both (comprehensive evaluation)
