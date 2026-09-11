#!/bin/bash
#SBATCH --job-name=run_load_balancing_1000_samples_5_clusters      # Job name
#SBATCH --output=./logs/%x-%j.out    # Standard output file
#SBATCH --error=./logs/%x-%j.err    # Standard error file
#SBATCH --partition=short     # Partition/queue name
#SBATCH --nodes=1              # Number of nodes/machines 
#SBATCH --ntasks=1             # Number of tasks/separate processes
#SBATCH --cpus-per-task=4     # CPU cores per task 
#SBATCH --mem=128G               # amount of ram 
#SBATCH --time=10:00:00        # Time limit hrs:min:sec

## run args ## 
export ix_sample_ids=500
export lb_sample_ids=1000
export threads=4
export n_clusters=5
export kmer_size=20
export base_expirements="/scratch/w.galbraith/expirements"
export scaled_num=1000

fastq_list="/scratch/w.galbraith/expirements/expirement_fastq_list/load_balance_files_from_${lb_sample_ids}_sample_ids.txt"
expirement_dir="/scratch/w.galbraith/expirements/${lb_sample_ids}_sample_ids_and_${n_clusters}_clusters"
intermediate_sketch_dir="/scratch/w.galbraith/expirements/1500_sample_ids_and_5_clusters/intermediate_sketches"
rr_sketch_dir="${expirement_dir}/round_robin_sketches"
wr_sketch_dir="${expirement_dir}/weighted_random_sketches"
sim_sketch_dir="${expirement_dir}/similarity_sketches"
initial_index_sketch_dir="${expirement_dir}/initial_index_sketches"

mkdir -p $rr_sketch_dir $sim_sketch_dir $wr_sketch_dir
cp $fastq_list $expirement_dir

## run the methods ## 

cargo run --release \
        -- run-similarity \
	-l $fastq_list \
	-i $intermediate_sketch_dir \
        -d $initial_index_sketch_dir \
        -t $threads \
        -o $sim_sketch_dir \
        -s $scaled_num \
        -k $kmer_size

# ## base lines ##
cargo run --release \
        -- run-round-robin \
	-l $fastq_list \
        -d $initial_index_sketch_dir \
        -o $rr_sketch_dir \
        -s $scaled_num \
        -k $kmer_size
cargo run --release \
        -- run-weighted-random \
        -o $wr_sketch_dir \
        -a "${sim_sketch_dir}/similarity_assignments.csv"
