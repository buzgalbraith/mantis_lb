#!/bin/bash
mkdir -p ./sra_data/fastq_files
curl https://www.cs.cmu.edu/~ckingsf/software/bloomtree/srr-list.txt | awk '{print $2}' > sra_data/accessions.txt
N=$1
tail -n "$N" sra_data/accessions.txt | while read acc; do
    echo "Fetching $acc..."
    prefetch "$acc" --output-directory sra_data/prefetch_files/
    fasterq-dump sra_data/prefetch_files/"$acc"/"$acc".sra \
        --outdir sra_data/fastq_files \
        --threads 8 \
        --split-files
done