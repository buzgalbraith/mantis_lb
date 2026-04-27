#!/bin/bash
#SBATCH --job-name=build_initial_sketch_clusters      # Job name
#SBATCH --output=./logs/%x-%j.out    # Standard output file
#SBATCH --error=./logs/%x-%j.err    # Standard error file
#SBATCH --partition=short     # Partition/queue name
#SBATCH --nodes=1              # Number of nodes/machines 
#SBATCH --ntasks=1             # Number of tasks/separate processes
#SBATCH --cpus-per-task=16      # CPU cores per task 
#SBATCH --mem=32G               # amount of ram 
#SBATCH --time=12:00:00        # Time limit hrs:min:sec

## run args ## 
export fastq_dir="/scratch/w.galbraith/CS7800_group_4/mantis/sra_data/initial_index_fastq_files"
export write_dir="/scratch/w.galbraith/CS7800_group_4/mantis/report_output/sketches/initial_index"
export n_clusters=5
export kmer_size=20
export scaled_num=1000

## run the method ## 
cargo run --release \
        -- build-index \
	-d $fastq_dir \
        -o $write_dir \
        -n $n_clusters \
        -s $scaled_num \
        -k $kmer_size
