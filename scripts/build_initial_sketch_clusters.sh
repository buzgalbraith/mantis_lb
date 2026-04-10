#!/bin/bash
#SBATCH --job-name=build_sketch_clusters      # Job name
#SBATCH --output=./logs/%x-%j.out    # Standard output file
#SBATCH --error=./logs/%x-%j.err    # Standard error file
#SBATCH --partition=short     # Partition/queue name
#SBATCH --nodes=1              # Number of nodes/machines 
#SBATCH --ntasks=1             # Number of tasks/separate processes
#SBATCH --cpus-per-task=16      # CPU cores per task 
#SBATCH --mem=32G               # amount of ram 
#SBATCH --time=02:00:00        # Time limit hrs:min:sec

## run the method ## 
cargo run --release \
        -- build-index \
        -d /scratch/w.galbraith/CS7800_group_4/mantis/sra_data/100_examples/ \
        -o ./initial_index \
        -n 5 \
        -s 1000 \
        -k 21
