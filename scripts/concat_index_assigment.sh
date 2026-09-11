#!/bin/bash
#SBATCH --job-name=concat_index_assignment
#SBATCH --output=./logs/%x-%A_%a.out    # %A = array job ID, %a = task index
#SBATCH --error=./logs/%x-%A_%a.err
#SBATCH --partition=short     # Partition/queue name
#SBATCH --nodes=1              # Number of nodes/machines
#SBATCH --ntasks=1             # Number of tasks/separate processes
#SBATCH --cpus-per-task=32      # CPU cores per task
#SBATCH --mem=64G               # amount of ram
#SBATCH --time=00:30:00        # Time limit hrs:min:sec

export kmer_size=20
export slots=31
export threads=32
export n_clusters=5
export base_expirements="/scratch/w.galbraith/expirements"
export scaled_num=1000
export base_squeakr="/scratch/w.galbraith/squeakr_files/kmer_size_20"
export num_ix_sample_ids=500 # number of sample ids for initial index
export num_lb_sample_ids=250 # number of sample ids for total index

## step 1 -> get assignment file for initial load balance files ## 
expirement_dir="${base_expirements}/${num_lb_sample_ids}_sample_ids_and_${n_clusters}_clusters"
ix_fastq_list="${base_expirements}/expirement_fastq_list/initial_index_files_from_${num_ix_sample_ids}_sample_ids.txt"
assignment_dir="${expirement_dir}/assignment_dir"
initial_index_sketch_dir="${expirement_dir}/initial_index_sketches"
mkdir -p $assignment_dir
cargo run --release \
        -- run-round-robin \
        -l $ix_fastq_list \
        -d $initial_index_sketch_dir \
        -o $assignment_dir \
        -s $scaled_num \
        -k $kmer_size

## step 2 -> concat initial index and load ballance assignments ## 
rr_assignment_file="${expirement_dir}/round_robin_sketches/round_robin_assignments.csv"
wr_assignment_file="${expirement_dir}/weighted_random_sketches/weighted_random_assignments.csv"
sim_assignment_file="${expirement_dir}/similarity_sketches/similarity_assignments.csv"
ix_assignment_file="${assignment_dir}/initial_index_assignment.csv"
mv ${assignment_dir}/round_robin_assignments.csv ${ix_assignment_file}
cp $rr_assignment_file ${assignment_dir}/round_robin_assignments.csv
cat ${ix_assignment_file} | grep '.fastq' >> ${assignment_dir}/round_robin_assignments.csv
cp $wr_assignment_file ${assignment_dir}/weighted_random_assignments.csv
cat ${ix_assignment_file} | grep '.fastq' >> ${assignment_dir}/weighted_random_assignments.csv
cp $sim_assignment_file ${assignment_dir}/similarity_assignments.csv
cat ${ix_assignment_file} | grep '.fastq' >> ${assignment_dir}/similarity_assignments.csv
## re-define for ease
rr_assignment_file="${assignment_dir}/round_robin_assignments.csv"
wr_assignment_file="${assignment_dir}/weighted_random_assignments.csv"
sim_assignment_file="${assignment_dir}/similarity_assignments.csv"
