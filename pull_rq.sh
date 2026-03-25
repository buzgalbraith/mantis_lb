#!/bin/zsh 
dir_name="/scratch/w.galbraith/CS7800_group_4/mantis/sra_data/temp"
mkdir -p fastq_files
for fname in SRR014494.fastq  SRR014495.fastq  SRR037700.fastq  SRR037701.fastq  SRR037702.fastq  SRR037703.fastq  SRR037705.fastq  SRR037706.fastq  SRR037707.fastq  SRR037708.fastq; do
	echo $fname
	rsync explorer:$dir_name/$fname fastq_files
done
