# Status Badge for README

Add this badge to the top of your `README.md` to show CI test status:

## Markdown Format

```markdown
[![CI Tests](https://github.com/HudoGriz/SV_coding_regions_benchmark_nextflow/actions/workflows/ci.yml/badge.svg)](https://github.com/HudoGriz/SV_coding_regions_benchmark_nextflow/actions/workflows/ci.yml)
```

## Preview

[![CI Tests](https://github.com/HudoGriz/SV_coding_regions_benchmark_nextflow/actions/workflows/ci.yml/badge.svg)](https://github.com/HudoGriz/SV_coding_regions_benchmark_nextflow/actions/workflows/ci.yml)

## Badge States

The badge will show:
- 🟢 **passing** - All tests successful
- 🔴 **failing** - One or more tests failed
- 🟡 **running** - Tests currently executing
- ⚪ **no status** - No workflow runs yet

## Example README.md

```markdown
# SV Calling and Benchmarking Pipeline

[![CI Tests](https://github.com/HudoGriz/SV_coding_regions_benchmark_nextflow/actions/workflows/ci.yml/badge.svg)](https://github.com/HudoGriz/SV_coding_regions_benchmark_nextflow/actions/workflows/ci.yml)

A Nextflow DSL2 pipeline for calling structural variants across multiple sequencing technologies and benchmarking results with Truvari.

[rest of your README content...]
```

## Additional Badges (Optional)

You can also add:

### Nextflow Version
```markdown
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-brightgreen.svg)](https://www.nextflow.io/)
```

### License (if applicable)
```markdown
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
```

### Last Commit
```markdown
[![GitHub last commit](https://img.shields.io/github/last-commit/HudoGriz/SV_coding_regions_benchmark_nextflow)](https://github.com/HudoGriz/SV_coding_regions_benchmark_nextflow/commits/master)
```

### Code Size
```markdown
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/HudoGriz/SV_coding_regions_benchmark_nextflow)](https://github.com/HudoGriz/SV_coding_regions_benchmark_nextflow)
```

## Full Badge Section Example

```markdown
# SV Calling and Benchmarking Pipeline

[![CI Tests](https://github.com/HudoGriz/SV_coding_regions_benchmark_nextflow/actions/workflows/ci.yml/badge.svg)](https://github.com/HudoGriz/SV_coding_regions_benchmark_nextflow/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Nextflow DSL2 pipeline for calling structural variants across multiple sequencing technologies and benchmarking results with Truvari.
```

This creates a nice professional look with multiple status indicators!
