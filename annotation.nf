nextflow.enable.dsl=2

process ANNOTATION {
    tag "SnpEff Annotation"
    input:
        tuple val(sample), path(vcf_in)
    output:
        tuple val(sample), path("svisionpro/${sample}_sv_results_annotated.vcf")
    script:
    """
    snpEff ann Lmajor $vcf_in > svisionpro/${sample}_sv_results_annotated.vcf
    """
}

workflow {
    Channel.fromFilePairs("svisionpro/*.vcf", checkIfExists: true)
        .ifEmpty { error "No VCF files found for annotation!" }
        .set { vcf_files }
    
    vcf_files
        | map { tuple(it[0].baseName, it[1]) }
        | ANNOTATION()
}

