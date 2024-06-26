#!/usr/bin/env python

import pysam
import pandas as pd
import edlib
import mappy

UP_seq = "TCTTCAGCGTTCCCGAGA"
UP_seq_revcomp = "TCTCGGGAACGCTGAAGA"

const_3prime_full_forward = "TCAGACGTGTGCTCTTCCGATCT"  # right
const_3prime_full_reverse = "AGATCGGAAGAGCACACGTCTGA"  # right
const_5prime_full_forward = "CTACACGACGCTCTTCCGATCT"  # left

print(const_3prime_full_forward)
right_const = const_3prime_full_reverse[:20]
left_const = const_5prime_full_forward[-20:]

N_read_extract = 100000


def extract_trimmed_fastq_pairs(ultima_fastq, R1_fastq, R2_fastq):
    i = 0
    max_dist = 3

    timmed_length = 55

    R1 = open(R1_fastq, "w")
    R2 = open(R2_fastq, "w")

    with pysam.FastxFile(ultima_fastq) as R:
        for r in R:
            i += 1

            seq = r.sequence
            rlen = len(seq)

            # reconstruction libraries are in this range
            if rlen > 135 and rlen < 170:
                # find TruSeqR1 position in the first 50nt of the read

                begin_seq = seq[:50]
                accept_r1 = False
                pos_con_in_begin = begin_seq.find(left_const)
                if pos_con_in_begin >= 0:
                    accept_r1 = True
                    pos_con_in_begin += len(left_const)
                else:
                    edit = edlib.align(
                        left_const, begin_seq, "HW", "locations", max_dist
                    )
                    dist = edit["editDistance"]
                    if dist >= 0:
                        accept_r1 = True
                        locs = edit["locations"][0]
                        pos_con_in_begin = locs[1] + 1

                # find TruSeqR2 position in the last 50nt of the read

                end_seq = seq[-50:]
                accept_r2 = False
                pos_con_in_end = end_seq.find(right_const)
                if pos_con_in_end >= 0:
                    accept_r2 = True
                    dist = 0
                else:
                    edit = edlib.align(
                        right_const, end_seq, "HW", "locations", max_dist
                    )
                    dist = edit["editDistance"]
                    if dist >= 0:
                        accept_r2 = True
                        locs = edit["locations"][0]
                        pos_con_in_end = locs[0]

                if accept_r2 and accept_r1:
                    qual = r.quality

                    # Actually trim and split to get R1
                    trim_begin = pos_con_in_begin
                    r1_seq = seq[trim_begin : trim_begin + timmed_length]
                    r1_qual = qual[trim_begin : trim_begin + timmed_length]

                    # this is a very specfic insertion pattern which I manually replace!!

                    if r1_seq[8:15] == "TCCTTCA":
                        r1_seq = r1_seq[:9] + r1_seq[10:]
                        r1_qual = r1_qual[:9] + r1_qual[10:]

                    # Actually trim and split to get R2

                    trim_end = 50 - pos_con_in_end
                    r2_seq = mappy.revcomp(seq[-trim_end - timmed_length : -trim_end])
                    r2_qual = qual[-trim_end - timmed_length : -trim_end][::-1]

                    R1.write(f"@{r.name}_1\n")
                    R1.write(f"{r1_seq}\n")
                    R1.write("+\n")
                    R1.write(f"{r1_qual}\n")

                    R2.write(f"@{r.name}_2\n")
                    R2.write(f"{r2_seq}\n")
                    R2.write("+\n")
                    R2.write(f"{r2_qual}\n")

            # if i > N_read_extract:
            #    break

    R1.close()
    R2.close()


def seq_counter(seq_dict, seq_instance):
    if seq_dict.get(seq_instance) is None:
        seq_dict[seq_instance] = 1
    else:
        seq_dict[seq_instance] += 1


def quad_dict_store(quad_dict, quad_key, quad_items):
    if quad_dict.get(quad_key) is None:
        quad_dict[quad_key] = [quad_items]
    else:
        quad_dict[quad_key].extend([quad_items])


def edit_match(input_seq, target_seq, max_dist):
    if input_seq == target_seq:
        dist = 0
        match = True
    else:
        edit = edlib.align(input_seq, target_seq, "NW", "path", max_dist)
        dist = edit["editDistance"]
        if dist >= 0 and dist <= max_dist:
            cigar = edit["cigar"]
            if "D" in cigar or "I" in cigar:
                match = False
                dist = "indel"
            else:
                match = True
        else:
            match = False

    return (match, dist)


def quality_calc(seq, quals, bases_dict, quals_dict):
    for i in range(len(seq)):
        if bases_dict.get(str(i)) is None:
            bases_dict[str(i)] = {}
            seq_counter(bases_dict[str(i)], seq[i])
        else:
            seq_counter(bases_dict[str(i)], seq[i])

        if quals_dict.get(str(i)) is None:
            quals_dict[str(i)] = {}
            seq_counter(quals_dict[str(i)], quals[i])
        else:
            seq_counter(quals_dict[str(i)], quals[i])


def quality_df(quals_dict):
    quals_df = pd.DataFrame(quals_dict)
    quals_df = quals_df.T
    quals_df = quals_df.fillna(0)
    quals_df = quals_df.stack()
    quals_df = quals_df.reset_index()
    quals_df.columns = ["position", "value", "ncount"]
    quals_df.position = quals_df.position.astype("int") + 1
    counts_df = quals_df.groupby("position").sum()
    quals_df["position_cnt"] = quals_df.position.apply(
        lambda x: counts_df.loc[x].ncount
    )
    quals_df["freq"] = quals_df.ncount / quals_df.position_cnt * 100
    return quals_df
