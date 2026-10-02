# Changelog

## 2.0.0 (2026-10-02)

The release analysed in the revised manuscript "How Fragmented Interval Subsetting Alters Structural Variant Benchmarking". With the `study` profile, it produces every value reported in the Results and the Supplementary Tables.

### Calling and input preparation
- **Delly 1.7.3:** added for Illumina WGS. It is also added for the GRCh37 Illumina WES data, restricted to the capture targets by adding every base outside them to Delly's exclusion list, the equivalent of Manta's exome call regions.
- **GRCh38 BAM preparation:** `preparation/build_grch38_analysis_bams.sh` now removes only alignments on contigs absent from the analysis reference, their mates, and SA-tag entries that point to those contigs. The previous flag filter dropped discordant pairs, mate-unmapped reads and supplementary alignments.
- **0/0 calls:** calls genotyped homozygous reference are removed before benchmarking (`--exclude_homref_calls`, on by default).
- **Study BEDs:** the gene-panel BEDs and the Agilent SureSelect Human All Exon V5 capture BED now ship in `data/`, so a fresh clone can prepare every target.

### Benchmarking and analyses
- **Sensitivity benchmarks** (`--sensitivity_benchmarks`):
  - matching-threshold grid (`--refdist`, `--pctsize`, `--pctseq`);
  - containment against any overlap;
  - target padding, and candidate-side `--extend`.
- **Post-hoc analyses** (`--posthoc_analyses`), run as pipeline tasks:
  - metrics collector and run-to-run comparison (`--posthoc_compare_results`);
  - SV-type accounting;
  - recall and precision decomposition with the false-positive origin audit;
  - stratification after matching;
  - block bootstrap and rarefied ranks;
  - composition standardisation;
  - simulation fidelity;
  - precision without inversions;
  - truth-membership counts;
  - SVanalyzer as a second comparator (bioconda image, `--svanalyzer_container`);
  - every value the manuscript and its supplementary tables report (`posthoc/manuscript/`).
- **Transition audit:** allele-aware record identity, so per-haplotype truth copies are kept apart.
- **`study` profile:** sets every option the study used.

### Running
- **One driver for every run:** `bin/run_benchmark.sh`. It records the pipeline file hashes so the code can be identified on nodes without git.
- **No site assumptions:** the repository holds no scheduler, site or hardware settings. Containers are published, signed images pulled on demand.

## 1.0.0

The version used for the submitted manuscript.
