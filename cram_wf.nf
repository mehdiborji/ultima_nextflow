
params.cram_file = "/Users/mborji/ultima/crams/sample_fail.cram"

process INDEX {

    cpus 10

    input:
    path input_cram

    script:
    """
    samtools index -@ ${task.cpus} $input_cram
    """
}

process CRAMTOFASTQ {

    cpus 10

    input:
    path input_cram
    val lines

    output:
    path 'fastq_chunk_*'

    script:
    """
    samtools fastq -@ ${task.cpus} $input_cram | split  -d -a 3 -l $lines - fastq_chunk_
    """
}

workflow {

    index_ch = INDEX(params.cram_file)
    fastq_chunks_ch = CRAMTOFASTQ(params.cram_file, 20000000)

}