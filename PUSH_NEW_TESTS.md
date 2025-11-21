# Instructions to Push PacBio and ONT Test Changes

## Summary of Changes

I've successfully added comprehensive testing for PacBio and ONT SV callers to your pipeline! Here's what was created:

### Files Added/Modified

1. **conf/test_pacbio.config** (NEW)
   - Test configuration for PacBio SV callers (CuteSV and PBSV)
   - Uses nf-core test datasets with PacBio HiFi data
   - Resource-limited for CI environments

2. **conf/test_ont.config** (NEW)
   - Test configuration for ONT SV callers (CuteSV and Sniffles)
   - Uses nf-core test datasets with Oxford Nanopore data
   - Includes tandem repeat annotations for Sniffles

3. **.github/workflows/ci.yml** (MODIFIED)
   - Added `test_pacbio` and `test_ont` to profile validation matrix
   - New dedicated test job: `test-pacbio` (tests CuteSV and PBSV)
   - New dedicated test job: `test-ont` (tests CuteSV and Sniffles)
   - Output validation and artifact collection for each technology

4. **nextflow.config** (MODIFIED)
   - Registered `test_pacbio` and `test_ont` profiles in profiles section

5. **docs/TESTING_LONG_READ_CALLERS.md** (NEW)
   - Comprehensive documentation for testing long-read SV callers
   - Usage examples and expected outputs
   - Troubleshooting guide
   - CI/CD integration details

### What Gets Tested

#### PacBio Test (`test_pacbio`)
- ✅ CuteSV on PacBio HiFi data
- ✅ PBSV discover + call workflow
- ✅ VCF compression and indexing
- ✅ Output file generation

#### ONT Test (`test_ont`)
- ✅ CuteSV on ONT data
- ✅ Sniffles with tandem repeat support
- ✅ VCF compression and indexing
- ✅ Output file generation

## How to Push Your Changes

Since this is based on PR #7 branch, you need to push to your fork and update the PR. Here are the steps:

### Option 1: Push to PR #7 Branch Directly (Recommended)

If you want to add these tests directly to PR #7:

```bash
cd SV_coding_regions_benchmark_nextflow

# Switch back to the PR #7 branch
git checkout pr-7

# Merge the new test branch
git merge add-pacbio-ont-tests

# Push to your fork
git push origin pr-7
```

This will automatically update PR #7 with the new tests.

### Option 2: Create a Separate PR for Tests

If you prefer a separate PR for the testing infrastructure:

```bash
cd SV_coding_regions_benchmark_nextflow

# Push the new branch to your fork
git push origin add-pacbio-ont-tests
```

Then go to GitHub and create a PR from `add-pacbio-ont-tests` to PR #7's target branch.

### Option 3: Create New PR to Main Repository

```bash
cd SV_coding_regions_benchmark_nextflow

# Push to your fork
git push origin add-pacbio-ont-tests
```

Then create a PR to the main repository.

## Testing Locally Before Pushing

You can test the new profiles locally to ensure everything works:

### Test PacBio Callers
```bash
nextflow run . -profile test_pacbio,docker --outdir test_results_pacbio
```

### Test ONT Callers
```bash
nextflow run . -profile test_ont,docker --outdir test_results_ont
```

### Validate Profile Loading
```bash
nextflow config -profile test_pacbio
nextflow config -profile test_ont
```

## What Happens After Pushing

Once you push and create/update a PR:

1. **GitHub Actions will automatically run** three test suites:
   - Existing Illumina tests (test_nfcore profile)
   - New PacBio tests (test_pacbio profile) - **YOUR NEW TESTS**
   - New ONT tests (test_ont profile) - **YOUR NEW TESTS**

2. **Profile validation** will check all 4 test profiles load correctly:
   - test
   - test_nfcore
   - test_pacbio ← NEW
   - test_ont ← NEW

3. **Artifacts will be collected** for each test run:
   - test-results (Illumina)
   - test-results-pacbio (PacBio) ← NEW
   - test-results-ont (ONT) ← NEW

## Expected CI Results

The CI will run these jobs:
- ✅ lint
- ✅ nf-core-lint
- ✅ validate-schema
- ✅ module-check
- ✅ profile-check (now tests 4 profiles)
- ✅ run-test (Illumina)
- ✅ test-pacbio (CuteSV + PBSV) ← NEW
- ✅ test-ont (CuteSV + Sniffles) ← NEW

## Troubleshooting

If CI fails, check:

1. **Test data availability**: Ensure nf-core test datasets are accessible
2. **Container images**: Docker/Singularity images need to be available
3. **Logs**: Download the workflow artifacts from GitHub Actions
4. **Local testing**: Run tests locally first with `-profile test_pacbio,docker`

## Documentation

For complete documentation on testing long-read callers, see:
- `docs/TESTING_LONG_READ_CALLERS.md`

## Questions?

If you have any questions about these changes or need modifications:
- Check the comprehensive docs in `docs/TESTING_LONG_READ_CALLERS.md`
- Review the CI workflow in `.github/workflows/ci.yml`
- Examine the test configs in `conf/test_pacbio.config` and `conf/test_ont.config`

## Changes Overview

```
5 files changed, 485 insertions(+), 1 deletion(-)

New files:
- conf/test_pacbio.config (68 lines)
- conf/test_ont.config (68 lines)
- docs/TESTING_LONG_READ_CALLERS.md (244 lines)

Modified files:
- .github/workflows/ci.yml (added ~100 lines)
- nextflow.config (added 8 lines)
```

## Commit Message

The commit includes a detailed message explaining:
- Test profiles added
- CI/CD enhancements
- Documentation created
- Features and benefits

Ready to push! 🚀
