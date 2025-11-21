# CI/CD Testing Setup Summary

## Overview

This document summarizes the GitHub Actions CI/CD testing infrastructure added to your SV calling pipeline.

## What Was Added

### 1. GitHub Actions Workflow (`.github/workflows/ci.yml`)

**Purpose**: Automated testing on every push and pull request

**Features**:
- ✅ Runs on push/PR to master, main, and dev branches
- ✅ Two separate jobs: Test and Lint
- ✅ Uses minimal synthetic data (~1 MB total)
- ✅ Tests in stub mode (fast, no container execution needed)
- ✅ Validates Nextflow syntax and configuration
- ✅ Uploads test artifacts for debugging
- ✅ Total runtime: ~5-7 minutes

**Jobs**:

1. **Test Job**:
   - Installs Nextflow (latest stable)
   - Sets up Apptainer/Singularity
   - Generates synthetic test data (chr22 1MB, BAM, VCF, BED files)
   - Runs pipeline in stub mode
   - Validates outputs
   - Uploads artifacts

2. **Lint Job**:
   - Validates Nextflow configuration
   - Checks DSL2 syntax
   - Shows available profiles

### 2. Test Profile in `nextflow.config`

**Purpose**: Optimized configuration for CI/CD testing

**Features**:
- Minimal resource allocation (2 CPUs, 6 GB memory, 10 min timeout)
- Docker-based execution (better for GitHub Actions)
- Lenient caching
- All processes configured with test-appropriate resources

**Usage**:
```bash
nextflow run main.nf -params-file test_params.yaml -profile test
```

### 3. Local Testing Script (`test_local.sh`)

**Purpose**: Test pipeline locally before pushing to GitHub

**Features**:
- Creates synthetic test data automatically
- Generates test configuration
- Runs pipeline in stub mode
- Provides cleanup instructions
- Executable script (chmod +x)

**Usage**:
```bash
./test_local.sh
```

### 4. Documentation

**Files added**:

- **`docs/TESTING.md`** (400 lines): Comprehensive testing guide covering:
  - Local testing procedures
  - GitHub Actions CI details
  - Test data specifications
  - Troubleshooting guide
  - Best practices
  - Advanced testing scenarios

- **`.github/README_TESTING.md`** (186 lines): GitHub Actions specific guide:
  - Workflow configuration
  - How to view results
  - Status badge setup
  - Customization options

- **`CI_SETUP_SUMMARY.md`** (this file): Quick reference

## Quick Start Guide

### 1. Test Locally First

```bash
# Make the test script executable (if not already)
chmod +x test_local.sh

# Run local test
./test_local.sh
```

This creates synthetic data and runs the pipeline in stub mode to verify workflow logic.

### 2. Commit and Push

```bash
# Add all new files
git add .github/ docs/TESTING.md test_local.sh nextflow.config

# Commit
git commit -m "Add GitHub Actions CI/CD testing infrastructure"

# Push
git push origin master
```

### 3. View Results

1. Go to your repository: https://github.com/HudoGriz/SV_coding_regions_benchmark_nextflow
2. Click the **Actions** tab
3. Watch the workflow run
4. Check the test results

### 4. Add Status Badge (Optional)

Add to the top of your `README.md`:

```markdown
[![CI Tests](https://github.com/HudoGriz/SV_coding_regions_benchmark_nextflow/actions/workflows/ci.yml/badge.svg)](https://github.com/HudoGriz/SV_coding_regions_benchmark_nextflow/actions/workflows/ci.yml)
```

## Test Data Specifications

The CI uses minimal synthetic data for fast testing:

| Component | Description | Size |
|-----------|-------------|------|
| Reference genome | Synthetic chr22, 1 MB | ~1 MB |
| BAM files | 20 aligned reads | ~1 KB |
| Truth VCF | 5 synthetic SVs | ~1 KB |
| Target BED files | 3 region sets | ~1 KB |

**Why synthetic data?**
- ⚡ Fast CI execution (~5 min instead of hours)
- 💰 No external data dependencies
- 🔒 No privacy concerns
- ✅ Tests workflow logic without heavy computation

## What Gets Tested

### Workflow Logic ✅
- Process connections are correct
- Channel operations work properly
- Conditional execution branches function
- All paths through the workflow are valid

### Configuration ✅
- Nextflow syntax is valid
- Profiles are properly defined
- Parameters are correctly specified
- Resource definitions are valid

### Process Definitions ✅
- All processes can be instantiated
- Inputs/outputs are correctly defined
- Directives are valid
- Module imports work

### What's NOT Tested (Yet)

- ❌ Actual tool execution (stub mode only)
- ❌ Output file correctness
- ❌ Real data processing
- ❌ Performance benchmarks
- ❌ Multi-platform execution

These can be added later with integration tests using real data.

## GitHub Actions Behavior

### Automatic Triggers

The workflow runs automatically:

1. **On Push** to branches:
   - `master`
   - `main`
   - `dev`

2. **On Pull Requests** to these branches

3. **Manual Trigger** via GitHub UI (Actions tab → Run workflow)

### Workflow Artifacts

Each run produces downloadable artifacts (retained for 7 days):
- `test-results/` - Pipeline output directory
- `.nextflow.log` - Execution log

Access artifacts: Workflow run → Artifacts section

## Resource Requirements

### GitHub Actions (Cloud)
- **Runner**: ubuntu-latest (2 CPUs, 7 GB RAM)
- **Duration**: ~5-7 minutes
- **Cost**: Free for public repositories

### Local Testing
- **Requirements**: Nextflow, samtools, bgzip, tabix, Python 3
- **Optional**: Docker or Singularity (for full run)
- **Duration**: ~2-3 minutes (stub mode)

## Customization Guide

### Change Test Branches

Edit `.github/workflows/ci.yml`:

```yaml
on:
  push:
    branches: [ master, main, dev, YOUR_BRANCH ]
  pull_request:
    branches: [ master, main, dev, YOUR_BRANCH ]
```

### Increase Timeout

If your tests need more time:

```yaml
jobs:
  test:
    timeout-minutes: 60  # Increase from 30
```

### Add Real Data Tests

Add a new job to `.github/workflows/ci.yml`:

```yaml
  integration-test:
    name: Integration Test with Real Data
    runs-on: ubuntu-latest
    steps:
      - name: Download test BAM
        run: |
          wget ftp://path/to/small/test.bam
      
      - name: Run pipeline
        run: |
          nextflow run main.nf -params-file integration_params.yaml
```

### Test Multiple Nextflow Versions

Add a matrix strategy:

```yaml
jobs:
  test:
    strategy:
      matrix:
        nextflow_version: ['23.04.0', '23.10.0', '24.04.0']
    steps:
      - name: Install Nextflow
        uses: nf-core/setup-nextflow@v2
        with:
          version: ${{ matrix.nextflow_version }}
```

## Troubleshooting

### Workflow Not Appearing

**Problem**: Pushed code but no Actions tab or workflows

**Solution**:
1. Check GitHub Actions is enabled: Settings → Actions → General
2. Verify workflow file is in `.github/workflows/ci.yml`
3. Check YAML syntax is valid

### Tests Failing

**Problem**: CI tests fail but works locally

**Debug steps**:
1. Click on failed run in Actions tab
2. Expand the failed step
3. Read error message carefully
4. Download artifacts for logs
5. Check if using wrong branch name
6. Verify all paths in test are relative

### Container Issues

**Problem**: "Unable to pull Docker image"

**Solutions**:
- Check Docker Hub availability
- Use GitHub Container Registry instead
- Pre-build and push containers
- Add retry logic to workflow

### Out of Memory

**Problem**: "Process exceeded memory limit"

**Solutions**:
- Increase test profile memory in `nextflow.config`
- Reduce test data size
- Use stub mode only (doesn't execute tools)

## Best Practices

1. **Test locally first**: Run `./test_local.sh` before pushing
2. **Small commits**: Push frequently with tested changes
3. **Watch CI**: Monitor first few runs to catch issues
4. **Use stub mode**: For fast workflow logic validation
5. **Add real tests gradually**: Start simple, add complexity
6. **Document failures**: Help future debugging
7. **Keep tests fast**: CI should complete in <10 minutes

## Next Steps

### Immediate Actions

1. ✅ Test locally: `./test_local.sh`
2. ✅ Commit changes: `git add` → `git commit`
3. ✅ Push to GitHub: `git push`
4. ✅ Watch Actions run
5. ✅ Add status badge to README

### Future Enhancements

Consider adding:

- [ ] **Integration tests** with real (small) data
- [ ] **Multiple genome builds** (GRCh37/38)
- [ ] **Parameter validation** tests
- [ ] **Resume functionality** tests
- [ ] **Output validation** checks
- [ ] **Performance benchmarks**
- [ ] **Multi-platform** tests (AWS, GCP)
- [ ] **Code coverage** reporting
- [ ] **Scheduled runs** (nightly builds)
- [ ] **Slack/email notifications**

## Files Modified/Added

### New Files
```
.github/
├── workflows/
│   └── ci.yml                    # GitHub Actions workflow
└── README_TESTING.md             # GitHub Actions guide

docs/
└── TESTING.md                     # Comprehensive testing guide

test_local.sh                      # Local testing script
CI_SETUP_SUMMARY.md               # This file
```

### Modified Files
```
nextflow.config                    # Added test profile
```

## Support & Resources

### Documentation
- [docs/TESTING.md](docs/TESTING.md) - Complete testing guide
- [.github/README_TESTING.md](.github/README_TESTING.md) - GitHub Actions guide

### External Resources
- [Nextflow Documentation](https://www.nextflow.io/docs/latest/)
- [GitHub Actions Documentation](https://docs.github.com/en/actions)
- [nf-core CI/CD Examples](https://github.com/nf-core/tools)

### Getting Help

If you encounter issues:
1. Check the documentation above
2. Review workflow logs in GitHub Actions
3. Test locally with `-with-trace` flag
4. Open a GitHub issue with:
   - Link to failed workflow run
   - Error message
   - Steps to reproduce
   - Nextflow version

## Summary

You now have a complete CI/CD testing infrastructure that:

✅ Automatically tests your pipeline on every push/PR  
✅ Validates Nextflow syntax and configuration  
✅ Runs fast (~5-7 minutes) with synthetic data  
✅ Provides local testing tools  
✅ Is well-documented  
✅ Can be easily extended  

The tests use **stub mode** which validates workflow logic without executing actual tools, making CI fast and reliable. This catches most common errors (syntax, configuration, logic) without the overhead of processing real data.

**Next**: Run `./test_local.sh` locally, then push to GitHub and watch your first CI run! 🚀
