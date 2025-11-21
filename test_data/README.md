# Test Data for Simulation Testing

This directory contains minimal test data for validating the simulation and statistics features of the SV Coding Regions Benchmark pipeline.

## Test Data Files

### gencode_test.gtf.gz
A minimal GTF file containing:
- **3 genes** (TEST1, TEST2, TEST3) across 2 chromosomes
- **3 transcripts** with exons, CDS, and UTR annotations
- Sufficient structure to test target region simulation

**Genes:**
- `TEST1` on chr1:100000-110000 (+ strand)
  - 3 exons with CDS and UTR regions
- `TEST2` on chr1:200000-210000 (- strand)
  - 3 exons with CDS and UTR regions  
- `TEST3` on chr2:50000-55000 (+ strand)
  - 2 exons with CDS and UTR regions

## Test Configuration

### Profile: test_simulation

The `test_simulation` profile is designed to quickly test simulation features:

```bash
nextflow run main.nf -profile test_simulation,docker
```

**Test Parameters:**
- `simulate_targets = true` - Enables simulation
- `num_simulations = 5` - Small number for fast testing
- `gather_statistics = true` - Enables statistics generation
- Single technology (Illumina WES) for speed

### Expected Behavior

When running the test_simulation profile, the pipeline should:

1. ✅ Run standard SV calling with Illumina WES data
2. ✅ Perform Truvari benchmarking on standard targets
3. ✅ Simulate 5 target region sets using the test GTF
4. ✅ Benchmark each simulated target set
5. ✅ Generate statistics and plots from all results

### Expected Outputs

```
test_results_simulation/
├── manta/                          # Manta SV calling results
│   └── illumina_wes/
├── benchmarking/                   # Benchmarking results
│   └── illumina_wes/
│       └── Manta/
│           ├── high_confidence/
│           ├── gene_panel/
│           ├── wes_utr/
│           └── simulated/          # Simulated target results
│               ├── simulation_1/
│               ├── simulation_2/
│               ├── simulation_3/
│               ├── simulation_4/
│               └── simulation_5/
├── simulated_targets/              # Generated BED files
│   ├── simulation_1.bed
│   ├── simulation_2.bed
│   ├── simulation_3.bed
│   ├── simulation_4.bed
│   └── simulation_5.bed
└── statistics/                     # Analysis outputs
    ├── plots/
    │   └── *.png
    ├── tables/
    │   └── *.csv
    └── summary_statistics.txt
```

## Running Tests

### Quick Test (Docker)
```bash
nextflow run main.nf -profile test_simulation,docker
```

### Quick Test (Singularity)
```bash
nextflow run main.nf -profile test_simulation,singularity
```

### Test Simulation Only
To test just the simulation module without running the full pipeline:
```bash
nextflow run main.nf \
  -profile test_simulation,docker \
  --skip_benchmarking
```

### Test Statistics Only
To test statistics generation on existing results:
```bash
nextflow run main.nf \
  -profile test_simulation,docker \
  --simulate_targets false \
  --gather_statistics true
```

## Validation Checks

After running the test, verify:

1. **Simulation Files Created**
   ```bash
   ls test_results_simulation/simulated_targets/simulation_*.bed
   ```
   Should show 5 BED files

2. **Simulation Benchmarking Completed**
   ```bash
   ls test_results_simulation/benchmarking/illumina_wes/Manta/simulated/
   ```
   Should show 5 simulation result directories

3. **Statistics Generated**
   ```bash
   ls test_results_simulation/statistics/plots/
   ls test_results_simulation/statistics/tables/
   ```
   Should contain plot and table files

4. **Summary Statistics**
   ```bash
   cat test_results_simulation/statistics/summary_statistics.txt
   ```
   Should contain aggregated statistics

## Test Duration

Expected runtime on modest hardware (2 CPUs, 6GB RAM):
- **Full test**: ~5-10 minutes
- **Simulation only**: ~1-2 minutes
- **Statistics only**: ~30 seconds

## Troubleshooting

### Common Issues

**Issue**: GTF file not found
```
Solution: Ensure gencode_test.gtf.gz exists in test_data/
```

**Issue**: Python/R dependencies missing
```
Solution: Use Docker or Singularity profiles which have all dependencies
```

**Issue**: Insufficient memory
```
Solution: Increase max_memory in test_simulation.config or run on larger machine
```

**Issue**: Simulation produces no output
```
Solution: Check that GTF file has valid gene annotations with exon/CDS/UTR features
```

## Integration with CI/CD

This test can be integrated into continuous integration:

```yaml
# .github/workflows/test_simulation.yml
name: Test Simulation Features
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Setup Nextflow
        uses: nf-core/setup-nextflow@v1
      - name: Run simulation test
        run: |
          nextflow run main.nf -profile test_simulation,docker
```

## Creating Additional Test Data

If you need to create test data for the base pipeline (not just simulation):

1. **Reference FASTA** (genome.fa)
   ```bash
   # Create minimal reference with chr1 and chr2
   # Should match chromosomes in GTF file
   ```

2. **Truth VCF** (truth.vcf.gz)
   ```bash
   # Create minimal VCF with a few test SVs
   ```

3. **Target BED files**
   ```bash
   # high_conf.bed, gene_panel.bed, wes_utr.bed
   # Should cover regions with test genes
   ```

4. **BAM file** (illumina_wes.bam)
   ```bash
   # Minimal alignment file for testing SV callers
   ```

See the main test.config and test_nfcore.config for examples of full test datasets.

## Notes

- This is a **minimal** test dataset for development and validation
- For production benchmarking, use real data from GIAB or similar
- Test data is designed to run quickly, not to produce biologically meaningful results
- All test coordinates are arbitrary and for validation purposes only
