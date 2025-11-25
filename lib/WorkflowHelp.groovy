/*
========================================================================================
    Helper functions for the pipeline
========================================================================================
*/

class WorkflowHelp {

    //
    // Print help message
    //
    public static String helpMessage() {
        return """
        =====================================================
        SV CALLING AND BENCHMARKING PIPELINE
        =====================================================
        
        Usage:
          nextflow run main.nf -profile <profile> [options]
        
        Required Arguments:
          --fasta                Reference genome FASTA file
        
        Input BAM Files (at least one required):
          --illumina_wes_bam     Illumina WES BAM file
          --illumina_wgs_bam     Illumina WGS BAM file
          --pacbio_bam           PacBio BAM file
          --ont_bam              Oxford Nanopore BAM file
        
        Benchmarking (optional):
          --benchmark_vcf        Truth VCF for benchmarking
          --skip_benchmarking    Skip Truvari benchmarking (default: false)
          --high_confidence_targets  BED file with high confidence regions
          --gene_panel_targets   BED file with gene panel regions
          --wes_utr_targets      BED file with WES UTR regions
        
        Optional Arguments:
          --outdir               Output directory (default: results)
          --run_name             Run name (default: benchmarking_run)
          --tandem_repeats       Tandem repeats BED file (for Sniffles)
          
        Simulation Options:
          --simulate_targets     Enable target region simulation (default: false)
          --num_simulations      Number of simulations to run (default: 100)
          --gencode_gtf          GENCODE GTF annotation file for simulation
          
        Analysis Options:
          --gather_statistics    Generate statistics and plots (default: false)
        
        Profiles:
          test_nfcore            Run with nf-core test data
          test                   Run with minimal test data
          docker                 Use Docker containers
          singularity            Use Singularity containers
        
        Examples:
          # Run with test data
          nextflow run main.nf -profile test_nfcore,docker
          
          # Run with custom data (no benchmarking)
          nextflow run main.nf -profile docker \\
            --fasta ref.fa \\
            --illumina_wes_bam sample.bam \\
            --skip_benchmarking
          
          # Run with benchmarking
          nextflow run main.nf -profile docker \\
            --fasta ref.fa \\
            --illumina_wes_bam sample.bam \\
            --benchmark_vcf truth.vcf.gz \\
            --high_confidence_targets regions.bed \\
            --gene_panel_targets panel.bed \\
            --wes_utr_targets wes_utr.bed
        =====================================================
        """.stripIndent()
    }

    //
    // Validate input parameters
    //
    public static void validateParameters(params, log) {
        def errors = []

        // Check if at least one BAM is provided
        if (!params.illumina_wes_bam && !params.illumina_wgs_bam && 
            !params.pacbio_bam && !params.ont_bam) {
            errors << "No input BAM files specified. Please provide at least one BAM file."
        }

        // Check reference FASTA
        if (!params.fasta) {
            errors << "Reference FASTA file (--fasta) is required."
        }

        // Check benchmarking parameters
        if (!params.skip_benchmarking && params.benchmark_vcf) {
            if (!params.high_confidence_targets || !params.gene_panel_targets || !params.wes_utr_targets) {
                log.warn """
                =====================================================
                WARNING: Benchmarking enabled but not all target BED files provided.
                
                Provide all three target files for complete benchmarking:
                  --high_confidence_targets <bed>
                  --gene_panel_targets <bed>
                  --wes_utr_targets <bed>
                
                Or use --skip_benchmarking to skip benchmarking.
                =====================================================
                """.stripIndent()
            }
        }

        // Check simulation parameters
        if (params.simulate_targets) {
            if (!params.gencode_gtf) {
                errors << "Simulation enabled (--simulate_targets) but --gencode_gtf not provided."
            }
            if (!params.benchmark_vcf) {
                errors << "Simulation enabled (--simulate_targets) but --benchmark_vcf not provided."
            }
        }

        // Report errors if any
        if (errors.size() > 0) {
            log.error """
            =====================================================
            PARAMETER VALIDATION ERRORS:
            
            ${errors.collect { "  ✗ ${it}" }.join('\n')}
            =====================================================
            """.stripIndent()
            System.exit(1)
        }
    }

    //
    // Log pipeline information
    //
    public static void logPipelineInfo(params, log) {
        def technologies = []
        if (params.illumina_wes_bam) technologies << "Illumina WES (Manta)"
        if (params.illumina_wgs_bam) technologies << "Illumina WGS (Manta)"
        if (params.pacbio_bam) {
            if (params.skip_pbsv) {
                technologies << "PacBio (CuteSV only)"
            } else {
                technologies << "PacBio (CuteSV, PBSV)"
            }
        }
        if (params.ont_bam) technologies << "ONT (CuteSV, Sniffles)"

        log.info """
        =====================================================
        SV CALLING AND BENCHMARKING PIPELINE
        =====================================================
        Run name      : ${params.run_name}
        Reference     : ${params.fasta}
        Output dir    : ${params.outdir}
        
        Technologies:
        ${technologies.collect { "  ✓ ${it}" }.join('\n')}
        
        Benchmarking  : ${params.skip_benchmarking ? 'Disabled' : 'Enabled'}
        Simulation    : ${params.simulate_targets ? 'Enabled' : 'Disabled'}
        Statistics    : ${params.gather_statistics ? 'Enabled' : 'Disabled'}
        =====================================================
        """.stripIndent()
    }
}
