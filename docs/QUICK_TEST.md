# Quick Test Reference Card

## 🚀 One-Line Tests

### Full Simulation Test
```bash
nextflow run main.nf -profile test_simulation,docker
```

### Statistics Only
```bash
nextflow run main.nf -profile test_simulation,docker --simulate_targets false
```

### More Simulations
```bash
nextflow run main.nf -profile test_simulation,docker --num_simulations 10
```

## ✅ Quick Validation

```bash
# Check outputs exist
ls test_results_simulation/simulated_targets/*.bed | wc -l  # Should be 5

# Verify statistics
ls test_results_simulation/statistics/plots/*.png | wc -l   # Should be > 0

# Check benchmarking
find test_results_simulation/benchmarking -name "summary.json" | wc -l
```

## 📊 Expected Runtime

- Full test: ~5-10 minutes
- Simulation only: ~1-2 minutes  
- Statistics only: ~30 seconds

## 🐛 Quick Troubleshooting

| Issue | Quick Fix |
|-------|-----------|
| GTF not found | `ls test_data/gencode_test.gtf.gz` |
| Out of memory | Add `--max_memory 12.GB` |
| No outputs | Check `--simulate_targets true` |
| R script fails | Use Docker/Singularity profile |

## 📁 Expected Outputs

```
test_results_simulation/
├── simulated_targets/     # 5 BED files
├── benchmarking/          # Truvari results
└── statistics/            # Plots & tables
```

## 🔗 More Info

- Full Guide: `docs/TESTING_SIMULATION.md`
- Test Data: `test_data/README.md`
- CI Workflow: `.github/workflows/test_simulation.yml`
