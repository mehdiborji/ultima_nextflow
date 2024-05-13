#!/usr/bin/env python

import argparse
import cram_utils


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_fastq", type=str)
parser.add_argument("-r1", "--read1_fastq", type=str)
parser.add_argument("-r2", "--read2_fastq", type=str)

args = parser.parse_args()
input_fastq = args.input_fastq
read1_fastq = args.read1_fastq
read2_fastq = args.read2_fastq

cram_utils.extract_trimmed_fastq_pairs(input_fastq, read1_fastq, read2_fastq)
