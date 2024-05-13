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

    cpus 4

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


process TRIM_EXTRACT_FASTQ_PAIRS_USING_BIN {

    input:
    path fastq_chunk

    output:
    path('read_*.fastq', arity: '2')

    script:
    """
    #!/usr/bin/env python

    import sys
    import os
    bin_path = os.path.abspath('../../../bin')
    print(bin_path)
    sys.path.append(bin_path)
    import cram_utils
    cram_utils.extract_trimmed_fastq_pairs("$fastq_chunk", "read_1.fastq", "read_2.fastq")

    """

}

process TRIM_EXTRACT_FASTQ_PAIRS_USING_BIN_COMMAND {

    input:
    path fastq_chunk

    output:
    path('read_*.fastq', arity: '2')

    script:
    """
    test.py -i $fastq_chunk -r1 read_1.fastq -r2 read_2.fastq
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

process ENV_TEST {

    output: stdout

    script:
    """
    #!/usr/bin/env python
    import os
    print(os.getcwd())
    import sys
    bin_path = os.path.abspath('../../../bin')

    print("old_way",bin_path)

    bin_path = os.path.abspath("$projectDir/bin")

    print("new_way",bin_path)

    #sys.path.append(bin_path)
    #import cram_utils
    """
    

}

workflow {

    //index_ch = INDEX(input_cram_ch)
    //fastq_chunks_ch = CRAMTOFASTQ(input_cram_ch, 500000)
    //read_pairs_ch = TRIM_EXTRACT_FASTQ_PAIRS(fastq_chunks_ch.flatten())
    //read_pairs_ch = TRIM_EXTRACT_FASTQ_PAIRS_USING_BIN(fastq_chunks_ch.flatten())
    //read_pairs_ch = TRIM_EXTRACT_FASTQ_PAIRS_USING_BIN_COMMAND(fastq_chunks_ch.flatten())
    //read_pairs_ch = TEST_USING_BIN("Hello")
    //read_pairs_ch.view { it }
    //outcounts = COUNT_READS(read_pairs_ch)
    //outcounts.view { it }

    cwd = ENV_TEST()
    cwd.view { it }

}