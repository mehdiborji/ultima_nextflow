
params.cram_file = "/Users/mborji/ultima/crams/sample_fail.cram"

process INDEX {

    cpus 4

    input:
    path input_cram

    script:
    """
    samtools index -@ ${task.cpus} $input_cram
    """
}

workflow {
    
    index_ch = INDEX(params.cram_file)
}