# WES Sequencing Targets Explained

## Overview

When running Manta on **Illumina WES (Whole Exome Sequencing)** data, you need to provide the **exome capture targets** BED file. This parameter was missing and has now been added.

## Parameter: `wes_sequencing_targets`

### What is it?
The BED file containing the genomic regions targeted by your exome capture kit.

### Why is it needed?
- Manta uses this as the `--callRegions` parameter
- Restricts SV calling to exome regions only
- Improves performance (faster runtime)
- Improves accuracy (fewer false positives in off-target regions)
- Required for proper WES analysis

### How is it different from other target files?

| Parameter | Purpose | When Used |
|-----------|---------|-----------|
| **`wes_sequencing_targets`** | Manta `--callRegions` for WES | SV calling on WES data |
| `wes_utr_targets` | WES + UTR regions for benchmarking | Truvari benchmarking |
| `gene_panel_targets` | Gene panel regions for benchmarking | Truvari benchmarking |
| `high_confidence_targets` | GIAB high confidence regions | Truvari benchmarking |

---

## Usage

### In params file:

```yaml
# For Illumina WES
illumina_wes_bam: '/path/to/wes.bam'

# Provide your exome capture kit BED file
wes_sequencing_targets: '/path/to/exome_capture_targets.bed'
```

### Can it be null?

**Yes**, it's optional:
- If `null` or not provided → Manta runs in WGS mode (scans entire genome)
- If provided → Manta runs in WES mode (scans only exome regions)

**Recommendation:** Always provide it for WES data for better results.

---

## Finding Your Exome Capture BED File

### Common Exome Kits:

1. **Agilent SureSelect**
   - Often named: `SureSelect_*.bed`
   - Example: `SureSelect_Human_All_Exon_V6_r2.bed`

2. **Illumina TruSeq**
   - Often named: `TruSeq_*.bed`
   - Example: `TruSeq_Exome_Targeted_Regions_v1.2.bed`

3. **Twist Bioscience**
   - Often named: `Twist_*.bed`
   - Example: `Twist_Exome_RefSeq_targets.bed`

4. **IDT xGen**
   - Often named: `xGen_*.bed`
   - Example: `xGen_Exome_Research_Panel_v2.bed`

### Where to find it:

1. **From your sequencing facility** - They should provide it with the data
2. **Kit manufacturer website** - Download from Agilent, Illumina, Twist, IDT, etc.
3. **In your project directory** - Check `data/targets/` or similar
4. **From UCSC Genome Browser** - If using standard kits

### Format requirements:

```
chr1    12345   67890   target_1
chr1    70000   90000   target_2
chr2    10000   20000   target_3
```

Must be:
- Tab-separated
- Sorted by chromosome and position
- Same reference build as your BAM (GRCh37/hg19 or GRCh38/hg38)
- Optional: bgzipped with tabix index (.bed.gz + .bed.gz.tbi)

---

## Example Configuration

### Minimal WES Configuration:

```yaml
# params_wes_minimal.yaml
run_name: 'wes_sv_calling'
outdir: './wes_results'

# Reference
fasta: '/data/reference/hs37d5.fa'

# WES data
illumina_wes_bam: '/data/bams/sample_wes.bam'

# Exome capture targets (CRITICAL!)
wes_sequencing_targets: '/data/targets/SureSelect_V6_exome.bed'

# Required target files (can be dummy for SV calling only)
high_confidence_targets: '/data/targets/dummy.bed'
gene_panel_targets: '/data/targets/dummy.bed'
wes_utr_targets: '/data/targets/dummy.bed'

# Skip benchmarking if no truth set
skip_benchmarking: true
```

### Full WES with Benchmarking:

```yaml
# params_wes_full.yaml
run_name: 'wes_benchmarking'
outdir: './wes_results'

# Reference
fasta: '/data/reference/hs37d5.fa'

# WES data
illumina_wes_bam: '/data/bams/HG002_wes.bam'

# Exome capture targets for Manta
wes_sequencing_targets: '/data/targets/SureSelect_V6_exome.bed'

# Benchmarking
benchmark_vcf: '/data/giab/HG002_SVs_Tier1_v0.6.vcf.gz'
skip_benchmarking: false

# Benchmarking target regions
high_confidence_targets: '/data/giab/HG002_SVs_Tier1_v0.6.bed'
gene_panel_targets: '/data/targets/gene_panel.bed'
wes_utr_targets: '/data/targets/exome_plus_utr.bed'
```

---

## What Happens Behind the Scenes

### With `wes_sequencing_targets` provided:

```bash
# Manta command line will include:
configManta.py \
  --bam sample.bam \
  --reference reference.fa \
  --callRegions exome_capture_targets.bed \  # ← This restricts calling
  --runDir manta
```

**Result:** Manta only calls SVs in exome regions

### Without `wes_sequencing_targets`:

```bash
# Manta command line:
configManta.py \
  --bam sample.bam \
  --reference reference.fa \
  --runDir manta
```

**Result:** Manta scans entire genome (slower, more off-target calls)

---

## Troubleshooting

### Issue: "File not found" error

```
Error: Cannot find file: exome_capture_targets.bed
```

**Solution:** 
1. Check the path in params file is correct
2. Use absolute path: `/full/path/to/file.bed`
3. Make sure file exists: `ls -lh /path/to/exome_capture_targets.bed`

---

### Issue: Don't have the exome capture BED file

**Solution A - Contact your sequencing facility:**
```
Ask them: "What exome capture kit was used? 
Can you provide the target BED file?"
```

**Solution B - Use standard kit BED:**
```bash
# Download from manufacturer
# Example for Agilent SureSelect V6:
wget https://earray.chem.agilent.com/suredesign/...
```

**Solution C - Set to null (not recommended):**
```yaml
# This will make Manta scan entire genome
wes_sequencing_targets: null
```

---

### Issue: BED file has wrong chromosome names

Your BED uses "chr1" but reference uses "1" (or vice versa):

**Solution - Fix chromosome names:**
```bash
# If BED has "chr" prefix but reference doesn't:
sed 's/^chr//' exome_targets.bed > exome_targets_fixed.bed

# If BED lacks "chr" prefix but reference has it:
sed 's/^/chr/' exome_targets.bed > exome_targets_fixed.bed
```

---

### Issue: BED file is for wrong genome build

Your BED is GRCh38 but BAM is GRCh37 (or vice versa):

**Solution - Use liftover:**
```bash
# Use UCSC liftOver tool
liftOver input.bed hg38ToHg19.over.chain.gz output.bed unMapped
```

---

## Best Practices

1. ✅ **Always provide for WES data** - Don't skip this parameter
2. ✅ **Match your capture kit** - Use the exact BED file for your kit
3. ✅ **Check genome build** - Ensure BED matches BAM reference
4. ✅ **Verify chromosome names** - "chr1" vs "1" must match
5. ✅ **Keep BED files organized** - Store in `data/targets/` directory
6. ✅ **Document your kit** - Add kit name to file or params comments

---

## Summary

### Required Files for WES Analysis:

```
✅ Reference genome (hs37d5.fa)
✅ WES BAM file + index
✅ Exome capture targets BED ← NEW! This is wes_sequencing_targets
✅ Three benchmarking BED files (can be dummy if skip_benchmarking=true)
```

### Complete Minimal WES Checklist:

- [ ] Reference FASTA + index
- [ ] WES BAM + BAI index
- [ ] **Exome capture targets BED** ← Most commonly missing!
- [ ] high_confidence_targets BED (can be dummy)
- [ ] gene_panel_targets BED (can be dummy)
- [ ] wes_utr_targets BED (can be dummy)

Once you have all these, your WES pipeline will run correctly! 🎉
