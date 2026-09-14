#!/bin/bash
#SBATCH --job-name=build_initial_sketch_clusters     # Job name
#SBATCH --output=./logs/%x-%j.out    # Standard output file
#SBATCH --error=./logs/%x-%j.err    # Standard error file
#SBATCH --partition=short     # Partition/queue name
#SBATCH --nodes=1              # Number of nodes/machines 
#SBATCH --ntasks=1             # Number of tasks/separate processes
#SBATCH --cpus-per-task=4     # CPU cores per task 
#SBATCH --mem=200G               # amount of ram 
#SBATCH --time=17:00:00        # Time limit hrs:min:sec

## run args ## 
# export fastq_dir="/scratch/w.galbraith/CS7800_group_4/mantis/sra_data/initial_index_fastq_files"
export ix_sample_ids=500 ## number of initial index sample ids
export lb_sample_ids=1500 ## number of load balancing sample ids
export fastq_list="/scratch/w.galbraith/expirements/diverse_fastq_list/initial_index_files_from_45_sample_ids.txt"
export threads=4
export n_clusters=5
export kmer_size=20
export scaled_num=1000

expirement_dir="/scratch/w.galbraith/expirements/diverse/${lb_sample_ids}_sample_ids_and_${n_clusters}_clusters/"
intermediate_sketch_dir="${expirement_dir}/intermediate_sketches"
write_dir="${expirement_dir}/initial_index_sketches"
mkdir -p $write_dir
cp $fastq_list $expirement_dir

## run the method ## 
cargo run --release \
        -- sketch-initial-index \
	-l $fastq_list \
	-i $intermediate_sketch_dir \
        -t $threads \
        -o $write_dir \
        -n $n_clusters \
        -s $scaled_num \
        -k $kmer_size
