nextflow.enable.dsl=2

process ALIGNMENT {
    tag "BWA and Minimap2 Alignment"
    input:
        tuple val(sample), path(reads1), path(reads2), path(ref)
    output:
        tuple val(sample), path("aligned_bwa_sorted.bam"), path("aligned_minimap2_sorted.bam")
    script:
    """
    # BWA Alignment
    bwa mem -t 8 $ref $reads1 $reads2 | samtools view -bS - | samtools sort -o ${sample}_bwa_sorted.bam
    samtools index ${sample}_bwa_sorted.bam
    
    # Minimap2 Alignment
    minimap2 -ax sr $ref $reads1 $reads2 | samtools view -bS - | samtools sort -o ${sample}_minimap2_sorted.bam
    samtools index ${sample}_minimap2_sorted.bam
    """
}

workflow {
    Channel.fromPath("data/*.fastq.gz", checkIfExists: true)
        .ifEmpty { error "No FASTQ files found!" }
        .set { fastq_files }
    
    fastq_files
        | map { it.toList().sort() }
        | map { tuple(it[0].baseName, it[0], it[1], "reference.fasta") }
        | ALIGNMENT()
}

