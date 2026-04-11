#!/bin/bash
#SBATCH --job-name=build_mantis_index      # Job name
#SBATCH --output=./logs/%x-%j.out    # Standard output file
#SBATCH --error=./logs/%x-%j.err    # Standard error file
#SBATCH --partition=short     # Partition/queue name
#SBATCH --nodes=1              # Number of nodes/machines
#SBATCH --ntasks=1             # Number of tasks/separate processes
#SBATCH --cpus-per-task=32      # CPU cores per task
#SBATCH --mem=64G               # amount of ram
#SBATCH --time=04:00:00        # Time limit hrs:min:sec

## load modules ##
module purge
module load Boost/1.88.0
## run vars ##
export n=5 ## number of clusters ##
export cluster_assignment_file="/scratch/w.galbraith/CS7800_group_4/mantis/sra_data/toy_example_files/toy_cluster_assignment.csv" ## csv file with cluster assignment ## 
export kmer_size=28
export slots=31
export write_dir="./mantis_output/"
export squeakr_dir="./mantis_output/squeakr_files"
export threads=32 ## make sure this matches the resource request#

## bins/lib paths ##
export MANTIS_BIN="/scratch/w.galbraith/CS7800_group_4/mantis/bin/mantis"
## run method ##
# ## get a list of squakr files from a dir, we may want to change this later ## 
for i in $(seq 0 $n); do
	echo "Building cluster ${i} index"
	## get list of squawker files ## 
	export index_dir="${write_dir}/cluster_${i}_index"
	mkdir -p $index_dir
	## get the files assigned to that cluster, replace the name with the corresponding squakr file and write that to a list##
	awk -F ',' -v i="$i" '$2 == i {print $1}' $cluster_assignment_file | xargs -I {} basename {} .fastq | xargs -I {} echo "${squeakr_dir}/{}.squeakr"  > "${index_dir}/cluster_${i}_squeakr_list.csv"
	## build mantis index ##
	echo "building mantis index..."
	$MANTIS_BIN build \
		-s $slots \
		-i "${index_dir}/cluster_${i}_squeakr_list.csv" \
		-o $index_dir
	echo "Mantis index built"
	## build mantis mst ##
	echo "Building mantis mst"
	$MANTIS_BIN mst \
		-p $index_dir \
		-t $threads \
		-k
	echo "Mantis MST done."
done
