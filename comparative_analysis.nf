nextflow.enable.dsl=2

// Process for comparing structural variant calls between Manta and SvisionPro
process COMPARATIVE_ANALYSIS {
    tag "Comparative Analysis of Manta and SvisionPro"
    input:
        tuple val(sample), path(vcf_manta), path(vcf_annotated) // Input includes sample name, VCF file from Manta, and annotated VCF from SvisionPro
    output:
        tuple val(sample), path("comparison/${sample}_comparison_results.csv") // Output is a CSV file summarizing differences
    script:
    """
    # Run comparative analysis using an R script
    Rscript comparative_analysis.R $vcf_manta $vcf_annotated comparison/${sample}_comparison_results.csv
    """
}

// Workflow to process and compare VCF files from both variant callers
workflow {
    Channel.fromFilePairs("manta/results/*.vcf", checkIfExists: true)
        .ifEmpty { error "No Manta VCF files found!" } // Ensure that input files exist
        .set { manta_vcf }
    
    Channel.fromFilePairs("svisionpro/*.vcf", checkIfExists: true)
        .ifEmpty { error "No SvisionPro VCF files found!" } // Ensure that input files exist
        .set { svisionpro_vcf }
    
    manta_vcf
        | join(svisionpro_vcf, by: 0) // Match VCF files for the same sample
        | map { tuple(it[0].baseName, it[1], it[2]) } // Prepare inputs for comparative analysis
        | COMPARATIVE_ANALYSIS() // Run the comparative analysis process
}

