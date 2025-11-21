# How to Add the GitHub Actions Workflow File

## Why This is Needed

The GitHub Actions workflow file (`.github/workflows/ci.yml`) couldn't be pushed automatically because it requires a GitHub token with the `workflow` scope. This is a security feature to prevent unauthorized modification of CI/CD workflows.

## Option 1: Add via GitHub Web UI (Easiest)

1. **Go to your repository**: https://github.com/HudoGriz/SV_coding_regions_benchmark_nextflow

2. **Navigate to the branch**:
   - Click on the branch dropdown (shows "master")
   - Select `seqera-ai/20251120-134739-add-github-actions-ci-testing`

3. **Create the workflow directory**:
   - Click "Add file" → "Create new file"
   - In the filename box, type: `.github/workflows/ci.yml`
   - GitHub will automatically create the nested directories

4. **Copy the workflow content**:
   - Open the file in your local clone or from the sandbox
   - Copy the entire contents of `.github/workflows/ci.yml` (shown below)
   - Paste into the GitHub editor

5. **Commit directly to the branch**:
   - Add commit message: "Add GitHub Actions workflow file"
   - Select "Commit directly to the `seqera-ai/20251120-134739-add-github-actions-ci-testing` branch"
   - Click "Commit new file"

## Option 2: Add via Git with Proper Token

If you have a personal access token with `workflow` scope:

```bash
# Fetch the branch
git fetch origin
git checkout seqera-ai/20251120-134739-add-github-actions-ci-testing

# Create the workflow file
mkdir -p .github/workflows
cat > .github/workflows/ci.yml << 'EOF'
[paste the workflow content here]
EOF

# Commit and push
git add .github/workflows/ci.yml
git commit -m "Add GitHub Actions workflow file"
git push origin seqera-ai/20251120-134739-add-github-actions-ci-testing
```

## Option 3: Add After Merging PR

You can also:
1. Merge the PR #6 as-is
2. Create a new branch from master
3. Add the workflow file
4. Create a small follow-up PR

## Workflow File Content

The complete content of `.github/workflows/ci.yml` is available in the sandbox or can be found in the PR branch. Here's a summary of what it does:

### Test Job
- Installs Nextflow and Apptainer
- Generates synthetic test data (1MB)
- Runs pipeline in stub mode
- Validates outputs
- Uploads artifacts

### Lint Job
- Validates Nextflow syntax
- Checks configuration
- Shows available profiles

## Verifying It Works

After adding the workflow file:

1. Go to the **Actions** tab in your repository
2. You should see workflows starting to appear
3. The workflow will run automatically on:
   - Push to master, main, or dev branches
   - Pull requests to these branches
   - Manual trigger (workflow_dispatch)

## Testing the Workflow

Once the workflow file is added to the PR branch:

1. Make a small commit to the branch (or just add the workflow file)
2. The workflow will run automatically
3. Check the Actions tab to see the results
4. Green checkmark = tests passed! ✅

## Troubleshooting

### Workflow Not Running

**Check**:
- File is at exactly `.github/workflows/ci.yml`
- File has valid YAML syntax
- GitHub Actions is enabled (Settings → Actions → General)

### Workflow Fails

**Common issues**:
- Syntax errors in YAML
- Missing dependencies
- Timeout issues (increase timeout in workflow file)

**Debug**:
- Click on failed run
- Expand failed step
- Read error message
- Download artifacts for logs

## Need Help?

- Check `docs/TESTING.md` for comprehensive testing guide
- Review `.github/README_TESTING.md` for GitHub Actions details
- See `CI_SETUP_SUMMARY.md` for quick reference

## After Adding the Workflow

Once the workflow is in place:

1. ✅ Merge PR #6
2. ✅ Watch first CI run in Actions tab
3. ✅ Add status badge to README (see `BADGE.md`)
4. ✅ Run `./test_local.sh` to test locally anytime

The workflow is configured to run on every push and PR, so you'll have automatic testing from that point forward! 🚀
