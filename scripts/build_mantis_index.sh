#!/bin/bash
#SBATCH --job-name=build_mantis_index_1500_sample_ids_7_clusters
#SBATCH --output=./logs/%x-%A_%a.out    # %A = array job ID, %a = task index
#SBATCH --error=./logs/%x-%A_%a.err
#SBATCH --partition=short     # Partition/queue name
#SBATCH --nodes=1              # Number of nodes/machines
#SBATCH --ntasks=1             # Number of tasks/separate processes
#SBATCH --cpus-per-task=32      # CPU cores per task
#SBATCH --mem=64G               # amount of ram
#SBATCH --time=12:00:00        # Time limit hrs:min:sec
#SBATCH --array=0-2

## load modules ##
module purge
module load Boost/1.88.0

export kmer_size=20
export slots=31
export threads=32
export n_clusters=9
export base_expirements="/scratch/w.galbraith/expirements"
export base_squeakr="/scratch/w.galbraith/squeakr_files"
export num_ix_sample_ids=500 # number of sample ids for initial index
export num_lb_sample_ids=1500 # number of sample ids for total index
expirement_dir="${base_expirements}/${num_lb_sample_ids}_sample_ids_and_${n_clusters}_clusters"
assignment_dir="${expirement_dir}/assignment_dir"
squeakr_dir="${base_squeakr}/kmer_size_${kmer_size}"
declare -A index_to_lb_method=(
    [0]="round_robin"
    [1]="weighted_random"
    [2]="similarity"
)
idx=$SLURM_ARRAY_TASK_ID
export MANTIS_BIN="/scratch/w.galbraith/CS7800_group_4/mantis/bin/mantis"
method="${index_to_lb_method[${idx}]}"
base_index_dir="${expirement_dir}/${method}_index"
assignment_file="${assignment_dir}/${method}_assignments.csv"
echo running $method
for i in $(seq 0 $(($n_clusters-1))); do
	echo "Building ${method} cluster ${i} index"
	export index_dir="${base_index_dir}/cluster_${i}_index"
	mkdir -p $index_dir
	awk -F ',' -v i="$i" '$2 == i {print $1}' $assignment_file \
		| xargs -I {} basename {} .fastq \
		| xargs -I {} echo "${squeakr_dir}/{}.squeakr" \
		> "${index_dir}/cluster_${i}_squeakr_list.csv"

	echo $index_dir
	echo "Building mantis index..."
	$MANTIS_BIN build \
		-s $slots \
		-i "${index_dir}/cluster_${i}_squeakr_list.csv" \
		-o $index_dir
	echo "Mantis index built"
	echo "Building mantis MST..."
	$MANTIS_BIN mst \
		-p $index_dir \
		-t $threads \
		-k
done
echo "Mantis MST done."
