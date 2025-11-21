# GitHub Actions CI Setup

This directory contains the GitHub Actions workflow configuration for automated testing of the SV calling pipeline.

## Workflow File

- **`.github/workflows/ci.yml`** - Main CI/CD workflow

## What Gets Tested

### On Every Push/PR:

1. **Syntax Validation**
   - Nextflow DSL2 syntax
   - Configuration validity
   - Profile definitions

2. **Workflow Logic**
   - Process connections
   - Channel operations
   - Stub execution (mock run)

3. **Code Quality**
   - Linting
   - Configuration checks

## Quick Start

### Enable GitHub Actions

GitHub Actions should be automatically enabled for your repository. If not:

1. Go to repository Settings
2. Navigate to Actions → General
3. Enable "Allow all actions and reusable workflows"

### View Test Results

1. Go to the **Actions** tab in your GitHub repository
2. Click on any workflow run to see details
3. Download artifacts to inspect test outputs

### Add Status Badge to README

Add this to the top of your `README.md`:

```markdown
[![CI Tests](https://github.com/HudoGriz/SV_coding_regions_benchmark_nextflow/actions/workflows/ci.yml/badge.svg)](https://github.com/HudoGriz/SV_coding_regions_benchmark_nextflow/actions/workflows/ci.yml)
```

This will display a badge showing the current test status.

## Workflow Triggers

The CI workflow runs:

- ✅ On push to `master`, `main`, or `dev` branches
- ✅ On pull requests to these branches  
- ✅ Manually via workflow_dispatch

## Test Duration

- **Lint job**: ~1-2 minutes
- **Test job**: ~3-5 minutes
- **Total**: ~5-7 minutes

## Test Data

Tests use synthetic minimal data:
- 1 MB synthetic reference genome (chr22)
- 20 aligned reads in BAM format
- 5 synthetic structural variants
- 3 target region BED files

**Total size**: ~1 MB (fast download and execution)

## Customization

### Modify Test Branches

Edit `.github/workflows/ci.yml`:

```yaml
on:
  push:
    branches: [ master, main, dev, YOUR_BRANCH ]
  pull_request:
    branches: [ master, main, dev, YOUR_BRANCH ]
```

### Adjust Timeout

If tests take longer:

```yaml
jobs:
  test:
    timeout-minutes: 30  # Increase this value
```

### Add More Tests

Add new steps to the `ci.yml` workflow:

```yaml
- name: Your Custom Test
  run: |
    # Your test commands here
    echo "Running custom test..."
```

## Artifacts

Each test run produces:
- `test-results/` - Pipeline outputs
- `.nextflow.log` - Nextflow execution log

**Retention**: 7 days (configurable in workflow)

## Troubleshooting

### Workflow Not Running

**Check:**
- GitHub Actions is enabled in repository settings
- Workflow file is in `.github/workflows/` directory
- Branch names match trigger conditions

### Tests Failing

**Debug steps:**
1. Click on failed workflow run
2. Expand failed job
3. Review error messages in logs
4. Download artifacts for detailed inspection
5. Test locally with `./test_local.sh`

### Container Pull Issues

If Docker pulls fail:
- Check Docker Hub status
- Consider using GitHub Container Registry
- Add retry logic to workflow

## Local Testing Before Push

Always test locally before pushing:

```bash
# Quick stub test
./test_local.sh

# Full test with containers
nextflow run main.nf -params-file test_params.yaml -profile test
```

This catches most issues before CI runs.

## Resources

- [GitHub Actions Documentation](https://docs.github.com/en/actions)
- [Nextflow CI/CD Guide](https://www.nextflow.io/docs/latest/sharing.html#continuous-integration)
- [nf-core CI Examples](https://github.com/nf-core/tools)

## Future Improvements

Potential enhancements:
- [ ] Add integration tests with real data
- [ ] Test multiple Nextflow versions
- [ ] Add code coverage reporting
- [ ] Set up scheduled runs (nightly)
- [ ] Add Slack/email notifications
- [ ] Test on multiple OS (Linux, macOS)
- [ ] Add security scanning

## Support

For issues with CI:
1. Check workflow logs
2. Test locally
3. Review [docs/TESTING.md](../docs/TESTING.md)
4. Open a GitHub issue with:
   - Link to failed workflow run
   - Error messages
   - Steps to reproduce
