# Instructions to Push Changes to GitHub

## Option 1: Push to a new branch and create PR

```bash
cd /home/user/SV_coding_regions_benchmark_nextflow

# Create and checkout a new branch
git checkout -b feature/replace-bedtools-with-nfcore

# Add all changes
git add -A

# Commit the changes
git commit -m "Replace local bedtools_intersect with nf-core module

- Updated all workflow files to use nf-core bedtools/intersect module
- Modified nextflow.config to configure nf-core module properly
- Updated process calls to match nf-core signature (added chrom_sizes parameter)
- Changed output references from .out.bed to .out.intersect
- Added BEDTOOLS_MIGRATION.md documentation

Benefits:
- Removes dependency on local Singularity containers
- Adds version tracking
- Follows nf-core best practices
- Easier maintenance and updates"

# Push to GitHub
git push -u origin feature/replace-bedtools-with-nfcore
```

After pushing, go to GitHub and you'll see a button to "Create Pull Request".

## Option 2: Push directly to master (if you have permissions)

```bash
cd /home/user/SV_coding_regions_benchmark_nextflow

# Add all changes
git add -A

# Commit
git commit -m "Replace local bedtools_intersect with nf-core module"

# Push to master
git push origin master
```

## Files Modified:
- ✅ nextflow.config
- ✅ workflows/preparation/prepare_data_complete_grch37.nf
- ✅ workflows/preparation/prepare_data_complete_grch38.nf
- ✅ workflows/prepare_data_grch37.nf
- ✅ workflows/prepare_data_grch38.nf

## New Files:
- ✅ BEDTOOLS_MIGRATION.md
- ✅ PUSH_INSTRUCTIONS.md (this file)

## PR Description Template:

```
## Summary
Replaces local bedtools_intersect module with nf-core bedtools/intersect module across all workflows.

## Changes
- Updated module includes from local to nf-core paths
- Modified process calls to match nf-core signature
- Updated output references from `.out.bed` to `.out.intersect`
- Added proper configuration in nextflow.config
- Added documentation in BEDTOOLS_MIGRATION.md

## Benefits
- ✅ No dependency on local Singularity containers
- ✅ Automatic container pulling from public registries
- ✅ Built-in version tracking
- ✅ Follows nf-core best practices
- ✅ Easier maintenance and updates

## Testing
Tested with nf-core modules already installed in the repository.

## Files Can Be Removed
- `modules/local/bedtools_intersect.nf` (no longer referenced)
```
