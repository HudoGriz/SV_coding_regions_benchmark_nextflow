# Pipeline Refactoring Summary

## Overview
The pipeline has been successfully refactored to follow nf-core best practices, improving code organization, maintainability, and reusability.

## New Modular Structure

### Sub-workflows Created

#### 1. **prepare_references.nf** (workflows/)
**Purpose**: Handles all reference file preparation and indexing

**Inputs**:
- Parameters (fasta, benchmark_vcf, target BED files, tandem_repeats)

**Outputs**:
- `fasta` - Reference genome
- `fasta_fai` - FAI index
- `benchmark_vcf` - Truth VCF (optional)
- `benchmark_vcf_tbi` - Truth VCF index (optional)
- `targets` - Channel with [name, bed_file] tuples
- `tandem_repeats` - Tandem repeats BED (optional)

**Features**:
- Automatic FAI index creation if missing
- Remote file support (HTTP/HTTPS/FTP)
- Conditional benchmark VCF handling
- Clean channel organization

---

#### 2. **sv_calling.nf** (workflows/)
**Purpose**: Orchestrates all SV calling across technologies

**Inputs**:
- `ch_fasta` - Reference FASTA
- `ch_fasta_fai` - Reference FAI index
- `ch_tandem_repeats` - Tandem repeats BED (optional)

**Outputs**:
- `vcfs` - Channel with [meta, vcf, tbi] tuples for all called SVs

**Features**:
- **Illumina WES**: Manta germline calling
- **Illumina WGS**: Manta germline calling
- **PacBio**: CuteSV + PBSV (optional)
- **ONT**: CuteSV + Sniffles
- Automatic VCF compression and indexing
- Technology-specific metadata propagation
- Conditional tool execution based on params

---

#### 3. **benchmarking.nf** (workflows/)
**Purpose**: Benchmarks SV calls against truth set using Truvari

**Inputs**:
- `ch_vcfs` - VCF files to benchmark
- `ch_benchmark_vcf` - Truth VCF
- `ch_benchmark_vcf_tbi` - Truth VCF index
- `ch_targets` - Target regions
- `ch_fasta` - Reference FASTA
- `ch_fasta_fai` - Reference FAI index

**Outputs**:
- `summary` - Truvari summary JSON files

**Features**:
- Automatic Truvari parameter selection (WES vs WGS)
- Per-target benchmarking
- Metadata-driven configuration
- Technology and tool tracking

---

### Main Workflow (main.nf)

The main workflow has been significantly simplified:

**Before**: ~350 lines with inline module calls
**After**: ~280 lines with clean sub-workflow calls

**Key Changes**:
1. **Cleaner imports**: Sub-workflows instead of individual modules
2. **Simplified logic**: 
   ```groovy
   // Old way (80+ lines)
   MANTA_WES(...)
   CUTESV_PACBIO(...)
   PBSV_DISCOVER(...)
   PBSV_CALL(...)
   CUTESV_ONT(...)
   SNIFFLES(...)
   // ... channel manipulation ...
   
   // New way (3 lines)
   SV_CALLING(ch_fasta, ch_fasta_fai, ch_tandem_repeats)
   ```

3. **Better separation of concerns**:
   - Data preparation → `PREPARE_REFERENCES()`
   - SV calling → `SV_CALLING()`
   - Benchmarking → `BENCHMARKING()`

4. **Maintained functionality**:
   - All existing features preserved
   - Same outputs and behavior
   - Backward compatible with existing params

---

## Benefits of Refactoring

### 1. **Improved Maintainability**
- Each sub-workflow has a single, clear purpose
- Changes to SV calling logic contained in one file
- Easier to debug and test individual components

### 2. **Better Code Reusability**
- Sub-workflows can be imported by other pipelines
- Modules organized by functionality
- Follows nf-core module structure

### 3. **Enhanced Readability**
- Main workflow now reads like a high-level flowchart
- Clear data flow between stages
- Reduced code duplication

### 4. **Easier Testing**
- Individual sub-workflows can be tested independently
- Smaller test scopes
- Faster debugging cycles

### 5. **Scalability**
- Easy to add new SV callers (just modify sv_calling.nf)
- Simple to extend benchmarking capabilities
- Can add new sub-workflows without cluttering main.nf

---

## File Structure

```
SV_coding_regions_benchmark_nextflow/
├── main.nf                          # Main workflow (simplified)
├── workflows/
│   ├── prepare_references.nf        # Reference preparation sub-workflow
│   ├── sv_calling.nf                # SV calling orchestration sub-workflow
│   ├── benchmarking.nf              # Truvari benchmarking sub-workflow
│   ├── simulate_and_benchmark.nf    # Simulation sub-workflow (existing)
│   ├── analysis_and_plots.nf        # Analysis sub-workflow (existing)
│   ├── prepare_giab_resources.nf    # GIAB prep (existing)
│   └── preparation/
│       ├── prepare_data_complete_grch37.nf
│       └── prepare_data_complete_grch38.nf
├── modules/
│   ├── nf-core/                     # nf-core modules (unchanged)
│   └── local/                       # Custom modules (unchanged)
└── lib/
    └── WorkflowHelp.groovy          # Helper functions (existing)
```

---

## Migration Notes

### For Pipeline Users
- **No changes required** to existing commands
- All parameters work exactly as before
- Same outputs and results
- Profiles unchanged

### For Pipeline Developers
- Edit `workflows/sv_calling.nf` to add new SV callers
- Edit `workflows/benchmarking.nf` to modify benchmarking logic
- Edit `workflows/prepare_references.nf` to change reference handling
- Main workflow (`main.nf`) remains clean and minimal

---

## Testing Checklist

✅ Test data preparation modes still work
✅ All SV calling technologies function correctly
✅ Benchmarking produces same results
✅ Simulation workflow integration
✅ Analysis and plotting
✅ Parameter validation
✅ Help message display
✅ Profile compatibility (docker, singularity, test)

---

## Future Improvements

Potential enhancements enabled by this refactoring:

1. **Add new SV callers** easily in `sv_calling.nf`:
   - Delly
   - Lumpy
   - GRIDSS
   - Dysgu

2. **Create additional benchmarking sub-workflows**:
   - Precision-recall curves
   - Multi-caller ensemble analysis
   - Stratified benchmarking

3. **Modular testing suite**:
   - Unit tests for each sub-workflow
   - Integration tests for full pipeline
   - Performance benchmarking

4. **Enhanced metadata propagation**:
   - Sample-level tracking
   - Run statistics collection
   - Automated reporting

---

## Conclusion

The refactoring successfully transforms a monolithic main workflow into a modular, maintainable pipeline following nf-core best practices. The code is now:
- **Cleaner**: Easier to read and understand
- **Modular**: Sub-workflows with clear responsibilities
- **Maintainable**: Changes isolated to specific components
- **Extensible**: Simple to add new features
- **Professional**: Follows community standards

All functionality has been preserved while significantly improving code quality.
