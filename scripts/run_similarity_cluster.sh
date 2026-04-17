#!/bin/bash
#SBATCH --job-name=run_similarity      # Job name
#SBATCH --output=./logs/%x-%j.out    # Standard output file
#SBATCH --error=./logs/%x-%j.err    # Standard error file
#SBATCH --partition=short     # Partition/queue name
#SBATCH --nodes=1              # Number of nodes/machines 
#SBATCH --ntasks=1             # Number of tasks/separate processes
#SBATCH --cpus-per-task=16      # CPU cores per task 
#SBATCH --mem=32G               # amount of ram 
#SBATCH --time=12:00:00        # Time limit hrs:min:sec

## run args ## 
export fastq_dir="/scratch/w.galbraith/CS7800_group_4/mantis/sra_data/load_ballance_fastq_files"
export initial_index_dir="/scratch/w.galbraith/CS7800_group_4/mantis/report_output/sketches/initial_index"
export write_dir="/scratch/w.galbraith/CS7800_group_4/mantis/report_output/sketches/similarity_index"
export n_clusters=5
export kmer_size=28
export scaled_num=1000

## run the method ## 
cargo run --release \
        -- run-similarity \
        -d $fastq_dir \
        -e $initial_index_dir \
        -o $write_dir \
        -s $scaled_num \
        -k $kmer_size
