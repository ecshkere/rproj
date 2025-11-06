#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.hisat_fa     = '/media/leon/Polina/Genomes/hg38.chromFa/hg38_FOR_HISAT.fa'
params.hisat_index  = '/media/leon/DISK2/icig/grch38_snp_tran/genome_snp_tran'
params.gtf          = '/media/leon/Polina/Genomes/Homo_sapiens.GRCh38.115.gtf'
params.metadata     = 'full_metadata.csv'
params.reads        = '/media/leon/DISK2/icig/done/fastqs/*{1,2}.fq.gz'
params.outdir       = '/media/leon/DISK2/icig/done/results'
params.threads      = 4
params.more_threads = 20


process FASTQC {
    publishDir "${params.outdir}/QC", mode: 'copy'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("*_fastqc.html"), path("*_fastqc.zip")

    script:
    """
    fastqc ${reads}
    """
}

process FASTP {
    cpus 4
    memory '8GB'
    maxForks 10
    publishDir "${params.outdir}/fastp", mode: 'copy'

    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("*.trimmed.fq.gz"), emit: reads
    path("*.json"), emit: json
    path("*.html"), emit: html

    script:
    """
    fastp \\
        -i ${reads[0]} -I ${reads[1]} \\
        -o ${sample}_1.trimmed.fq.gz \\
        -O ${sample}_2.trimmed.fq.gz \\
        --detect_adapter_for_pe --trim_poly_g \\
        --thread ${task.cpus} \\
        --length_required 36 \\
        -j ${sample}.fastp.json -h ${sample}.fastp.html
    """
}

process HISAT2 {
    cpus 40
    memory '70GB'
    maxForks 1   
    publishDir "${params.outdir}/alignments", mode: 'copy'

    input:
    tuple val(sample), path(reads), val(assay_type)

    output:
    tuple val(sample), path("${sample}.bam"), val(assay_type), emit: bam
    tuple val(sample), path("${sample}.bam.bai"), val(assay_type), emit: bai

    script:
    """
    splicing_mode="--no-spliced-alignment" 
    if [ ${assay_type} == "RNASEQ" ]; then
        splicing_mode=""
    fi

    hisat2 -p ${task.cpus} -x ${params.hisat_index} \\
        -1 ${reads[0]} -2 ${reads[1]} \${splicing_mode} \\
        | samtools view -@ ${task.cpus} -b -q 30 -u - \\
        | samtools sort -@ ${task.cpus} -o ${sample}.bam

    samtools index -@ ${task.cpus} ${sample}.bam
    """
}

process FEATURE_COUNTS {
    cpus 20
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path(alignments)

    output:
    tuple path("counts.tsv"), path("counts.tsv.summary")

    script:
    """
    featureCounts -a ${params.gtf} -o counts.tsv -T ${task.cpus} -p -B ${alignments} -Q 30 --countReadPairs
    """
}

process GET_STATS {
    cpus 1
    memory '6GB'
    maxForks 14

    publishDir "${params.outdir}/stats", mode: 'copy'

    input:
    tuple val(sample), path(bam_file), val(assay_type)

    output:
    tuple val(sample), path("*.stat")

    script:
    """
    stats.pl ${sample} ${bam_file} ${params.hisat_fa}
    """
}

process MERGE_BAMS {
    cpus 1
    memory '1GB'
    maxForks 20
    tag "$patient"
    publishDir "${params.outdir}/merged", mode: 'copy'

    input:
    tuple val(patient), path(bams), val(assay_type)

    output:
    tuple val(patient), path("${assay_type}_${patient}.bam"), emit: bam
    path("${assay_type}_${patient}.bam.bai"), emit: bai

    script:
    """
    samtools merge ${assay_type}_${patient}.merged.bam ${bams}
    samtools merge -r ${assay_type}_${patient}.merged.bam ${assay_type}_${patient}.bam ## joining samples (?)
    rm ${assay_type}_${patient}.merged.bam
    samtools index ${assay_type}_${patient}.bam
    """
}

workflow {
    fastqs = Channel.fromFilePairs(params.reads)
    trimmed_fastqs = FASTP(fastqs).reads
    FASTQC(trimmed_fastqs)

    metadata = Channel.fromPath(params.metadata)
        .flatMap { file -> 
            file.readLines().drop(1).collect { line ->
                def (patient, assay_type, condition, sample) = line.split(',')
                tuple(sample, assay_type)
            }
        }
    reads_dct = trimmed_fastqs.join(metadata)
    alignments = HISAT2(reads_dct).bam

    GET_STATS(alignments)

    alignments.branch {
        rna: it[2] == 'RNASEQ'
        chip: it[2] == 'CHIPSEQ'
    }
    .set { assay }

    rna_alignments = assay.rna.map { sample, bam_path, assay_type -> bam_path }.collect()
    // chip_alignments = assay.chip.map { sample, bam_path, assay_type -> bam_path }.collect()
    
    FEATURE_COUNTS(rna_alignments)
}
