# Migration: SV Callers to nf-core Modules

## Overview

This document describes the comprehensive migration of all SV caller modules (PBSV, CuteSV, Manta, and Sniffles) from local implementations to their nf-core equivalents. This is the most significant migration in the standardization initiative, affecting the core SV calling functionality of the pipeline.

## Migration Summary

| Tool | Local Module | nf-core Module | Status |
|------|-------------|----------------|--------|
| **Manta** | `modules/local/manta.nf` | `modules/nf-core/manta/germline/main.nf` | ✅ Migrated |
| **CuteSV** | `modules/local/cutesv.nf` | `modules/nf-core/cutesv/main.nf` | ✅ Migrated |
| **PBSV** | `modules/local/pbsv.nf` | `modules/nf-core/pbsv/discover` + `call` | ✅ Migrated |
| **Sniffles** | `modules/local/sniffles.nf` | `modules/nf-core/sniffles/main.nf` | ✅ Migrated |

### New Helper Module

Created `modules/local/bgzip_tabix.nf` to compress and index VCF outputs from CuteSV and PBSV (which output uncompressed VCFs in nf-core versions).

---

## Technical Changes

### 1. MANTA Migration

#### Key Differences

| Aspect | Local Module | nf-core Module |
|--------|-------------|----------------|
| **Process Name** | `MANTA` | `MANTA_GERMLINE` |
| **Input** | `tuple [meta, bam, bai]`, `path fasta`, `path fai` | `tuple [meta, bam, bai, target_bed, target_bed_tbi]`, `tuple [meta2, fasta]`, `tuple [meta3, fai]`, `path config` |
| **Output** | `diploidSV.vcf.gz` + `.tbi` in results/variants/ | Multiple VCFs: diploid_sv, candidate_sv, candidate_small_indels |
| **Container** | Local singularity | `community.wave.seqera.io/library/manta_python` |
| **Working Dir** | `manta_run/` | `manta/` |

#### Code Changes

**Before (Local)**:
```groovy
include { MANTA as MANTA_WES } from './modules/local/manta'

MANTA_WES(
    ch_illumina_wes_bam,  // [meta, bam, bai]
    ch_fasta,              // path
    ch_fasta_fai           // path
)
```

**After (nf-core)**:
```groovy
include { MANTA_GERMLINE as MANTA_WES } from './modules/nf-core/manta/germline/main'

ch_illumina_wes_bam = Channel.value([
    [id: 'Illumina_WES', technology: 'Illumina_WES', tool: 'Manta'],
    file(params.illumina_wes_bam),
    file("${params.illumina_wes_bam}.bai"),
    [],  // target_bed (empty for WGS/WES without specific targets)
    []   // target_bed_tbi
])

MANTA_WES(
    ch_illumina_wes_bam,
    ch_fasta.map { f -> [[id: 'fasta'], f] },
    ch_fasta_fai.map { f -> [[id: 'fai'], f] },
    []  // config file (optional)
)

// Use diploid_sv output
ch_manta_vcf = MANTA_WES.out.diploid_sv_vcf.join(MANTA_WES.out.diploid_sv_vcf_tbi)
```

---

### 2. CUTESV Migration

#### Key Differences

| Aspect | Local Module | nf-core Module |
|--------|-------------|----------------|
| **Input** | `tuple [meta, bam, bai]`, `path fasta`, `path fai` | `tuple [meta, bam, bai]`, `tuple [meta2, fasta]` |
| **Output** | `HG002_hs37d5.vcf.gz` + `.tbi` | `*.vcf` (uncompressed) |
| **Container** | Local singularity (cutesv:2.1.1) | `biocontainers/cutesv:2.0.2` |
| **Compression** | Built-in (bgzip + tabix) | Requires separate BGZIP_TABIX step |
| **Work Dir** | `work_dir/` | Current directory (`.`) |

#### Code Changes

**Before (Local)**:
```groovy
include { CUTESV as CUTESV_PACBIO } from './modules/local/cutesv'

CUTESV_PACBIO(
    ch_pacbio_bam,  // [meta, bam, bai]
    ch_fasta,        // path
    ch_fasta_fai     // path
)
// Output: HG002_hs37d5.vcf.gz + .tbi
```

**After (nf-core + BGZIP_TABIX)**:
```groovy
include { CUTESV as CUTESV_PACBIO } from './modules/nf-core/cutesv/main'
include { BGZIP_TABIX as BGZIP_TABIX_CUTESV_PACBIO } from './modules/local/bgzip_tabix'

CUTESV_PACBIO(
    ch_pacbio_bam.map { meta, bam, bai -> 
        [[id: meta.id, technology: meta.technology, tool: 'CuteSV'], bam, bai]
    },
    ch_fasta.map { f -> [[id: 'fasta'], f] }
)

// Compress and index output
BGZIP_TABIX_CUTESV_PACBIO(CUTESV_PACBIO.out.vcf)
// Output: [meta, *.vcf.gz, *.vcf.gz.tbi]
```

---

### 3. PBSV Migration

#### Key Differences

| Aspect | Local Module | nf-core Module |
|--------|-------------|----------------|
| **Architecture** | Single process | **Two-step**: PBSV_DISCOVER → PBSV_CALL |
| **Input** | `tuple [meta, bam, bai]`, `path fasta`, `path fai` | Discover: `tuple [meta, bam]`, `tuple [meta2, fasta]`<br>Call: `tuple [meta, svsig]`, `tuple [meta2, fasta]` |
| **Output** | `HG002_hs37d5.vcf.gz` + `.tbi` + `HG002.svsig.gz` | Discover: `*.svsig.gz`<br>Call: `*.vcf` (uncompressed) |
| **Container** | Local singularity (pbsv:2.10.0) | `community.wave.seqera.io/library/pbsv:2.11.0` |
| **Compression** | Built-in | Requires BGZIP_TABIX |

#### Code Changes

**Before (Local - single process)**:
```groovy
include { PBSV } from './modules/local/pbsv'

PBSV(
    ch_pacbio_bam,  // [meta, bam, bai]
    ch_fasta,        // path
    ch_fasta_fai     // path
)
// Output: HG002_hs37d5.vcf.gz + .tbi, HG002.svsig.gz
```

**After (nf-core - two processes + compression)**:
```groovy
include { PBSV_DISCOVER } from './modules/nf-core/pbsv/discover/main'
include { PBSV_CALL } from './modules/nf-core/pbsv/call/main'
include { BGZIP_TABIX as BGZIP_TABIX_PBSV } from './modules/local/bgzip_tabix'

// Step 1: Discover SV signatures
PBSV_DISCOVER(
    ch_pacbio_bam.map { meta, bam, bai -> 
        [[id: meta.id, technology: meta.technology, tool: 'Pbsv'], bam]
    },
    ch_fasta.map { f -> [[id: 'fasta'], f] }
)

// Step 2: Call SVs from signatures
PBSV_CALL(
    PBSV_DISCOVER.out.svsig,
    ch_fasta.map { f -> [[id: 'fasta'], f] }
)

// Step 3: Compress and index
BGZIP_TABIX_PBSV(PBSV_CALL.out.vcf)
// Output: [meta, *.vcf.gz, *.vcf.gz.tbi]
```

---

### 4. SNIFFLES Migration

#### Key Differences

| Aspect | Local Module | nf-core Module |
|--------|-------------|----------------|
| **Input** | `tuple [meta, bam, bai]`, `path fasta`, `path fai`, `path tandem_repeats` | `tuple [meta, bam, bai]`, `tuple [meta2, fasta]`, `tuple [meta3, tandem_file]`, `val vcf_output`, `val snf_output` |
| **Output** | `HG002_hs37d5.vcf.gz` + `.tbi` | `*.vcf.gz` + `.tbi` (when vcf_output=true) |
| **Container** | Local singularity (sniffles:2.4) | `biocontainers/sniffles:2.4` |
| **Compression** | Built-in | Built-in (controlled by vcf_output flag) |
| **Tandem Repeats** | Optional with NO_FILE check | Tuple with empty channel if not provided |

#### Code Changes

**Before (Local)**:
```groovy
include { SNIFFLES } from './modules/local/sniffles'

SNIFFLES(
    ch_ont_bam,          // [meta, bam, bai]
    ch_fasta,            // path
    ch_fasta_fai,        // path
    ch_tandem_repeats    // path (or NO_FILE)
)
// Output: HG002_hs37d5.vcf.gz + .tbi
```

**After (nf-core)**:
```groovy
include { SNIFFLES } from './modules/nf-core/sniffles/main'

SNIFFLES(
    ch_ont_bam.map { meta, bam, bai -> 
        [[id: meta.id, technology: meta.technology, tool: 'Sniffles'], bam, bai]
    },
    ch_fasta.map { f -> [[id: 'fasta'], f] },
    params.tandem_repeats ? 
        ch_tandem_repeats.map { f -> [[id: 'tandem_repeats'], f] } :
        Channel.value([[id: 'null'], []]),
    true,   // vcf_output (enable VCF output with compression)
    false   // snf_output (disable SNF format output)
)

// Use vcf and tbi outputs
ch_sniffles_vcf = SNIFFLES.out.vcf.join(SNIFFLES.out.tbi)
```

---

## Helper Module: BGZIP_TABIX

### Purpose

Compresses uncompressed VCF files with bgzip and indexes them with tabix. Required for CuteSV and PBSV nf-core modules which output uncompressed VCFs.

### Implementation

```groovy
process BGZIP_TABIX {
    tag "${meta.id}"
    label 'process_low'
    
    container 'community.wave.seqera.io/library/htslib:1.21--ff8e28a189fbecaa'
    
    input:
    tuple val(meta), path(vcf)
    
    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.tool}"
    """
    bgzip -c ${vcf} > ${prefix}.vcf.gz
    tabix -p vcf ${prefix}.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bgzip: \$(bgzip --version 2>&1 | head -1 | sed 's/^.*bgzip //' | sed 's/ .*\$//')
        tabix: \$(tabix --version 2>&1 | head -1 | sed 's/^.*tabix //' | sed 's/ .*\$//')
    END_VERSIONS
    """
}
```

### Usage Pattern

```groovy
// For any nf-core module outputting uncompressed VCF
TOOL(inputs...)
BGZIP_TABIX(TOOL.out.vcf)
// Now you have: [meta, vcf.gz, vcf.gz.tbi]
```

---

## Configuration Updates

### Resource Allocation

```groovy
// Manta
withName: 'MANTA_GERMLINE' {
    cpus = { check_max(8, 'cpus') }
    memory = { check_max(32.GB, 'memory') }
    time = { check_max(12.h, 'time') }
}

// CuteSV
withName: 'CUTESV_.*' {
    cpus = { check_max(params.max_cpus, 'cpus') }
    memory = { check_max(64.GB, 'memory') }
    time = { check_max(8.h, 'time') }
}

// PBSV
withName: 'PBSV_DISCOVER' {
    cpus = 1
    memory = { check_max(8.GB, 'memory') }
    time = { check_max(4.h, 'time') }
}

withName: 'PBSV_CALL' {
    cpus = { check_max(params.max_cpus, 'cpus') }
    memory = { check_max(32.GB, 'memory') }
    time = { check_max(8.h, 'time') }
}

// Sniffles
withName: 'SNIFFLES' {
    cpus = { check_max(params.max_cpus, 'cpus') }
    memory = { check_max(32.GB, 'memory') }
    time = { check_max(8.h, 'time') }
}

// BGZIP_TABIX helper
withName: 'BGZIP_TABIX.*' {
    cpus = 2
    memory = 4.GB
    time = 2.h
}
```

### Important: No Container Overrides

❌ **Removed** container overrides that were using local singularity images  
✅ **Now** nf-core modules use their default public containers

This ensures:
- Reproducibility across environments
- No dependency on local container builds
- Automatic version tracking in versions.yml

---

## Channel Structure Changes

### VCF Output Collection

**Before**:
```groovy
ch_all_vcfs = Channel.empty()
if (params.illumina_wes_bam) {
    ch_all_vcfs = ch_all_vcfs.mix(MANTA_WES.out.vcf)
}
if (params.pacbio_bam) {
    ch_all_vcfs = ch_all_vcfs.mix(
        CUTESV_PACBIO.out.vcf,
        PBSV.out.vcf
    )
}
```

**After**:
```groovy
ch_all_vcfs = Channel.empty()
if (params.illumina_wes_bam) {
    // MANTA outputs separate vcf and tbi channels - join them
    ch_all_vcfs = ch_all_vcfs.mix(
        MANTA_WES.out.diploid_sv_vcf.join(MANTA_WES.out.diploid_sv_vcf_tbi)
    )
}
if (params.pacbio_bam) {
    // BGZIP_TABIX outputs tuple [meta, vcf, tbi]
    ch_all_vcfs = ch_all_vcfs.mix(
        BGZIP_TABIX_CUTESV_PACBIO.out.vcf,
        BGZIP_TABIX_PBSV.out.vcf
    )
}
if (params.ont_bam) {
    // SNIFFLES outputs vcf and tbi separately - join them
    ch_all_vcfs = ch_all_vcfs.mix(
        BGZIP_TABIX_CUTESV_ONT.out.vcf,
        SNIFFLES.out.vcf.join(SNIFFLES.out.tbi)
    )
}
```

---

## Benefits of Migration

### 1. No Local Container Dependencies ✅

**Before**: Required pre-built local singularity images
- `file:///singularity_images/manta_latest.sif`
- `file:///singularity_images/cutesv_latest.sif`
- `file:///singularity_images/pbsv_latest.sif`
- `file:///singularity_images/sniffles_latest.sif`

**After**: Uses public containers automatically
- Manta: `community.wave.seqera.io/library/manta_python`
- CuteSV: `biocontainers/cutesv:2.0.2`
- PBSV: `community.wave.seqera.io/library/pbsv:2.11.0`
- Sniffles: `biocontainers/sniffles:2.4`

### 2. Automatic Version Tracking ✅

All nf-core modules generate `versions.yml` files capturing:
- Tool versions
- Container versions
- Execution environment

### 3. Standardized Interfaces ✅

- Consistent meta map structure
- Predictable input/output patterns
- Better documentation and examples

### 4. Easy Updates ✅

```bash
# Update individual modules
nf-core modules update manta/germline
nf-core modules update cutesv
nf-core modules update pbsv/discover pbsv/call
nf-core modules update sniffles

# Or update all at once
nf-core modules update --all
```

### 5. Community Maintenance ✅

- Benefit from nf-core community bug fixes
- Security patches
- Performance improvements
- New features

### 6. Better Testing ✅

nf-core modules include:
- Stub runs for testing
- Example test data
- CI/CD integration

---

## Testing

### Test Individual SV Callers

```bash
# Test Manta on Illumina WES
nextflow run main.nf \
    --fasta /path/to/reference.fa \
    --illumina_wes_bam /path/to/wes.bam \
    --benchmark_vcf /path/to/truth.vcf.gz \
    --high_confidence_targets /path/to/targets.bed \
    --gene_panel_targets /path/to/panel.bed \
    --wes_utr_targets /path/to/wes_utr.bed

# Test CuteSV and PBSV on PacBio
nextflow run main.nf \
    --fasta /path/to/reference.fa \
    --pacbio_bam /path/to/pacbio.bam \
    --benchmark_vcf /path/to/truth.vcf.gz \
    --high_confidence_targets /path/to/targets.bed \
    --gene_panel_targets /path/to/panel.bed \
    --wes_utr_targets /path/to/wes_utr.bed

# Test CuteSV and Sniffles on ONT
nextflow run main.nf \
    --fasta /path/to/reference.fa \
    --ont_bam /path/to/ont.bam \
    --tandem_repeats /path/to/tandem_repeats.bed \
    --benchmark_vcf /path/to/truth.vcf.gz \
    --high_confidence_targets /path/to/targets.bed \
    --gene_panel_targets /path/to/panel.bed \
    --wes_utr_targets /path/to/wes_utr.bed
```

### Expected Output Structure

```
outdir/
├── calls/
│   ├── Illumina_WES/
│   │   └── sv/
│   │       └── manta/
│   │           ├── Illumina_WES.diploid_sv.vcf.gz
│   │           └── Illumina_WES.diploid_sv.vcf.gz.tbi
│   ├── PacBio/
│   │   └── sv/
│   │       ├── cutesv/
│   │       │   ├── PacBio_CuteSV.vcf.gz
│   │       │   └── PacBio_CuteSV.vcf.gz.tbi
│   │       └── pbsv/
│   │           ├── PacBio_Pbsv.vcf.gz
│   │           └── PacBio_Pbsv.vcf.gz.tbi
│   └── ONT/
│       └── sv/
│           ├── cutesv/
│           │   ├── ONT_CuteSV.vcf.gz
│           │   └── ONT_CuteSV.vcf.gz.tbi
│           └── sniffles/
│               ├── ONT_Sniffles.vcf.gz
│               └── ONT_Sniffles.vcf.gz.tbi
└── benchmarking/
    └── truvari/
        └── ...
```

---

## Files That Can Be Removed

After this PR is merged, the following local modules are no longer used:

1. ✅ `modules/local/manta.nf` - Replaced by `nf-core/manta/germline`
2. ✅ `modules/local/cutesv.nf` - Replaced by `nf-core/cutesv`
3. ✅ `modules/local/pbsv.nf` - Replaced by `nf-core/pbsv/discover` + `call`
4. ✅ `modules/local/sniffles.nf` - Replaced by `nf-core/sniffles`

### Remaining Local Modules (Appropriate)

These modules should remain local as they contain pipeline-specific logic:
- `modules/local/bgzip_tabix.nf` - **NEW** helper module for VCF compression
- `modules/local/truvari.nf` - Custom benchmarking wrapper
- `modules/local/create_target_beds.nf` - Custom BED generation
- `modules/local/download_*.nf` - Custom download utilities

---

## Troubleshooting

### Issue: "No such process 'MANTA'"

**Cause**: Process renamed from `MANTA` to `MANTA_GERMLINE`

**Solution**: Update all references:
```groovy
// Old
MANTA(inputs...)

// New
MANTA_GERMLINE(inputs...)
```

### Issue: "PBSV output not found"

**Cause**: PBSV is now two-step process

**Solution**: Ensure both steps are completed:
```groovy
PBSV_DISCOVER(bam, fasta)
PBSV_CALL(PBSV_DISCOVER.out.svsig, fasta)
BGZIP_TABIX_PBSV(PBSV_CALL.out.vcf)
```

### Issue: "CuteSV VCF not compressed"

**Cause**: nf-core CUTESV outputs uncompressed VCF

**Solution**: Add BGZIP_TABIX step:
```groovy
CUTESV(bam, fasta)
BGZIP_TABIX(CUTESV.out.vcf)
// Use BGZIP_TABIX.out.vcf
```

### Issue: "Sniffles not producing output"

**Cause**: vcf_output flag not set to true

**Solution**: Enable VCF output:
```groovy
SNIFFLES(
    bam,
    fasta,
    tandem_repeats,
    true,   // vcf_output - MUST be true
    false   // snf_output
)
```

### Issue: "Container not found"

**Cause**: Singularity may need to pull container

**Solution**: Pre-pull containers or ensure internet access:
```bash
# Pre-pull containers
nextflow run main.nf --help  # Triggers container pulls
```

---

## Migration Progress

- ✅ **bedtools** - Completed (PR #1)
- ✅ **gunzip** - Completed (PR #2)
- ✅ **samtools_faidx** - Completed (PR #3)
- ✅ **manta** - Completed (this PR)
- ✅ **cutesv** - Completed (this PR)
- ✅ **pbsv** - Completed (this PR)
- ✅ **sniffles** - Completed (this PR)

### Final Local Modules Status

✅ **Appropriate to keep**:
- `bgzip_tabix.nf` - Helper module for VCF compression
- `truvari.nf` - Custom benchmarking wrapper
- `create_target_beds.nf` - Custom BED generation
- `download_*.nf` - Custom download utilities

---

## Container Versions

| Tool | nf-core Version | Local Version (old) |
|------|-----------------|---------------------|
| **Manta** | configManta.py (from container) | dceoy/manta:latest |
| **CuteSV** | 2.0.2 | 2.1.1 |
| **PBSV** | 2.11.0 | 2.10.0 |
| **Sniffles** | 2.4 | 2.4 |

**Note**: CuteSV nf-core version is slightly older (2.0.2 vs 2.1.1). If specific features from 2.1.1 are required, you can override the container in `nextflow.config`:

```groovy
withName: 'CUTESV_.*' {
    container = 'biocontainers/cutesv:2.1.1--pyhdfd78af_0'
}
```

---

## Related Documentation

- [bedtools Migration](BEDTOOLS_MIGRATION.md)
- [gunzip Migration](GUNZIP_MIGRATION.md)
- [samtools_faidx Migration](SAMTOOLS_TABIX_MIGRATION.md)
- [nf-core modules documentation](https://nf-co.re/modules)

---

## Questions?

If you encounter issues or have questions:
1. Check the nf-core module documentation for specific modules
2. Review the example workflows in this repository
3. Verify channel structures match expected inputs
4. Open an issue on GitHub with details

**Migration completed**: 2025-11-20  
**Author**: Seqera AI
