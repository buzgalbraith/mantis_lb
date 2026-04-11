#!/bin/bash
#SBATCH --job-name=run_round_robin      # Job name
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
        -- run-round-robin \
        -d /scratch/w.galbraith/CS7800_group_4/mantis/sra_data/toy_example_files \
        -e ./initial_index \
        -o ./sketches/round_robin_sketches \
        -s 1000 \
        -k 21
