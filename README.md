# ultima_nextflow

## Some methods to process the cram files from Ultima Genomics for downstream jobs

- First we subsample the inial file to %1 and also include only reads in >=120 and <=180 length range:

```
samtools view -@12 -f 4 -s .01 034339-H8_3_rec-UGAv3-115-CTCGCGCAATGCGAT.cram | awk 'length($10) >= 120 && length($10) <= 180' | samtools view -@6 -b -o output_sub1p.cram
```
- Nextflow can be run by giving this cram file as input:

```
nextflow run cram_wf.nf --cram_file /Users/mborji/nf_ultima/recons/output_sub1p.cram
```
