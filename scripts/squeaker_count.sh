#!/bin/bash
#SBATCH --job-name=squeaker_run      # Job name
#SBATCH --output=./outputs/%x-%j.out    # Standard output file
#SBATCH --error=./outputs/%x-%j.err    # Standard error file
#SBATCH --partition=short     # Partition/queue name
#SBATCH --nodes=1              # Number of nodes/machines
#SBATCH --ntasks=1             # Number of tasks/separate processes
#SBATCH --cpus-per-task=16      # CPU cores per task
#SBATCH --mem=32G               # amount of ram
#SBATCH --time=02:00:00        # Time limit hrs:min:sec

## run the method ##
module purge 
module load Boost/1.88.0
FASTQ_DIR=/scratch/w.galbraith/CS7800_group_4/mantis/sra_data/load_ballance_fastq_files
WRITE_DIR=./squeaker_outputs
SQUEAKER_BIN=/scratch/w.galbraith/CS7800_group_4/mantis/bin/squeakr 
LD_LIBRARY_PATH=/scratch/w.galbraith/CS7800_group_4/mantis/lib/boost-1.84.0/
mkdir -p $WRITE_DIR
for fastq_file in $FASTQ_DIR/*.fastq; do
	base=$(basename "$fastq_file" .fastq)
	write_name=$WRITE_DIR/$base.squeakr
	echo $write_name
	$SQUEAKER_BIN count -e \
		--no-counts \
		--cutoff 3 \
		-k 31 \
		-s 33 \
		-t 16 \
		-o $write_name \
		$fastq_file
done
