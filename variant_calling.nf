nextflow.enable.dsl=2

process VARIANT_CALLING {
    tag "Manta and SvisionPro Variant Calling"
    input:
        tuple val(sample), path(bam_bwa), path(bam_minimap), path(ref)
    output:
        tuple val(sample), path("manta/results/${sample}_variants.vcf"), path("svisionpro/${sample}_sv_results.vcf")
    script:
    """
    # Run Manta
    mkdir -p manta/${sample}
    configManta.py --bam $bam_bwa --referenceFasta $ref --runDir manta/${sample}
    manta/${sample}/runWorkflow.py -m local -j 8
    
    # Run SvisionPro
    mkdir -p svisionpro
    svisionpro --bam $bam_bwa --reference $ref --output svisionpro/${sample}_sv_results.vcf --min_sv_size 50 --max_sv_size 500000 --min_qual 30 --threads 8
    """
}

workflow {
    Channel.fromFilePairs("data/*.bam", checkIfExists: true)
        .ifEmpty { error "No BAM files found!" }
        .set { bam_files }
    
    bam_files
        | map { tuple(it[0].baseName, it[0], it[1], "reference.fasta") }
        | VARIANT_CALLING()
}

