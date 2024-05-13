#!/usr/bin/env nextflow

params.cram_file = "/Users/mborji/ultima/crams/sample_fail.cram"
input_cram_ch = Channel.of(params.cram_file)

process INDEX {

    cpus 4

    input:
    path input_cram

    script:
    """
    samtools index -@ ${task.cpus} $input_cram
    """
}

process CRAMTOFASTQ {

    cpus 8

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

process EXTRACT_TRIMMED_FASTQ_PAIRS {

    input:
    path fastq_chunk

    output:
    path('read_*.fastq', arity: '2')

    script:
    """
    #!/usr/bin/env python
    import sys
    import os

    bin_path = os.path.abspath("$projectDir/bin")
    sys.path.append(bin_path)

    import cram_utils
    cram_utils.extract_trimmed_fastq_pairs("$fastq_chunk", "read_1.fastq", "read_2.fastq")
    """

}

process COUNT_READS {

    input:
    path('read_*.fastq', arity: '2')

    output: stdout

    script:
    """
    echo "Counting read1"
    samtools view -c read_1.fastq
    echo "Counting read2"
    samtools view -c read_2.fastq
    """

}

workflow {

    //index_ch = INDEX(input_cram_ch)
    fastq_chunks_ch = CRAMTOFASTQ(input_cram_ch, 4000000)
    read_pairs_ch = EXTRACT_TRIMMED_FASTQ_PAIRS(fastq_chunks_ch.flatten())
    outcounts = COUNT_READS(read_pairs_ch)
    outcounts.view { it }

}