#!/bin/bash

# Script to push bedtools nf-core migration changes to GitHub

set -e

echo "=================================================="
echo "Pushing Bedtools nf-core Migration Changes"
echo "=================================================="
echo ""

cd /home/user/SV_coding_regions_benchmark_nextflow

# Show current status
echo "📋 Current git status:"
git status --short
echo ""

# Ask user preference
echo "How would you like to push these changes?"
echo "1) Create a new branch and push (recommended for PR)"
echo "2) Push directly to current branch"
read -p "Enter your choice (1 or 2): " choice

if [ "$choice" == "1" ]; then
    # Create new branch
    BRANCH_NAME="feature/replace-bedtools-with-nfcore"
    echo ""
    echo "🌿 Creating new branch: $BRANCH_NAME"
    git checkout -b $BRANCH_NAME
    
    # Add and commit
    echo "➕ Adding all changes..."
    git add -A
    
    echo "💾 Committing changes..."
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
    
    echo "🚀 Pushing to GitHub..."
    git push -u origin $BRANCH_NAME
    
    echo ""
    echo "✅ Success! Branch pushed to GitHub."
    echo ""
    echo "Next steps:"
    echo "1. Go to: https://github.com/HudoGriz/SV_coding_regions_benchmark_nextflow"
    echo "2. Click on 'Compare & pull request' button"
    echo "3. Review and create the PR"
    
elif [ "$choice" == "2" ]; then
    # Push to current branch
    echo ""
    echo "➕ Adding all changes..."
    git add -A
    
    echo "💾 Committing changes..."
    git commit -m "Replace local bedtools_intersect with nf-core module"
    
    CURRENT_BRANCH=$(git branch --show-current)
    echo "🚀 Pushing to branch: $CURRENT_BRANCH"
    git push origin $CURRENT_BRANCH
    
    echo ""
    echo "✅ Success! Changes pushed to GitHub."
else
    echo "❌ Invalid choice. Exiting."
    exit 1
fi

echo ""
echo "=================================================="
echo "📁 Files modified:"
echo "  - nextflow.config"
echo "  - workflows/preparation/prepare_data_complete_grch37.nf"
echo "  - workflows/preparation/prepare_data_complete_grch38.nf"
echo "  - workflows/prepare_data_grch37.nf"
echo "  - workflows/prepare_data_grch38.nf"
echo ""
echo "📄 Files added:"
echo "  - BEDTOOLS_MIGRATION.md"
echo "  - PUSH_INSTRUCTIONS.md"
echo "  - push_changes.sh"
echo "=================================================="
