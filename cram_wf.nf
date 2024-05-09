#!/usr/bin/env nextflow

params.cram_file = "/Users/mborji/ultima/crams/sample_fail.cram"
input_cram_ch = Channel.of(params.cram_file)

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

process TRIM_EXTRACT_FASTQ_PAIRS {

    input:
    path fastq_chunk

    output:
    path R1_fastq

    script:
    """
    #!/usr/bin/env python
    import pysam
    i = 0

    N_read_extract = 100000

    timmed_length = 55

    R1 = open(R1_fastq, "w")

    with pysam.FastxFile(fastq_chunk) as R:
        for r in R:
            i += 1

            seq = r.sequence
            rlen = len(seq)

            # reconstruction libraries are in this range
            if rlen == 155:

                trim_begin = 20
                qual = r.quality
                r1_seq = seq[trim_begin : trim_begin + timmed_length]
                r1_qual = qual[trim_begin : trim_begin + timmed_length]

                R1.write(f"@{r.name}_1\n")
                R1.write(f"{r1_seq}\n")
                R1.write("+\n")
                R1.write(f"{r1_qual}\n")
            
            if i > N_read_extract:
                break

    R1.close()
    
    """
}

workflow {

    index_ch = INDEX(input_cram_ch)
    fastq_chunks_ch = CRAMTOFASTQ(input_cram_ch, 20000000)
    read_pairs_ch = TRIM_EXTRACT_FASTQ_PAIRS(fastq_chunks_ch.flatten())

}