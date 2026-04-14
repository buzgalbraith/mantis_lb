#!/bin/bash
#SBATCH --job-name=build_mantis_index_round_robin
#SBATCH --output=./logs/%x-%A_%a.out    # %A = array job ID, %a = task index
#SBATCH --error=./logs/%x-%A_%a.err
#SBATCH --partition=short     # Partition/queue name
#SBATCH --nodes=1              # Number of nodes/machines
#SBATCH --ntasks=1             # Number of tasks/separate processes
#SBATCH --cpus-per-task=1      # CPU cores per task
#SBATCH --mem=64G               # amount of ram
#SBATCH --time=12:00:00        # Time limit hrs:min:sec
#SBATCH --array=0-3

## load modules ##
module purge
module load Boost/1.88.0

## run vars ##
export cluster_assignment_file="/scratch/pathak.shr/round_robin_assignments.csv"
export kmer_size=28
export slots=31
export write_dir="/scratch/w.galbraith/CS7800_group_4/mantis/mantis_output/round_robin_index"
export squeakr_dir="/scratch/w.galbraith/CS7800_group_4/mantis/mantis_output/squeakr_files/"
export threads=32

## bins/lib paths ##
export MANTIS_BIN="/scratch/w.galbraith/CS7800_group_4/mantis/bin/mantis"

## use the array task ID as the cluster index ##
i=$SLURM_ARRAY_TASK_ID

echo "Building cluster ${i} index"

export index_dir="${write_dir}/cluster_${i}_index"
mkdir -p $index_dir

awk -F ',' -v i="$i" '$2 == i {print $1}' $cluster_assignment_file \
	| xargs -I {} basename {} .fastq \
	| xargs -I {} echo "${squeakr_dir}/{}.squeakr" \
	> "${index_dir}/cluster_${i}_squeakr_list.csv"

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
echo "Mantis MST done."
