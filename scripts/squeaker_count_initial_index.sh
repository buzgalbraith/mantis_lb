#!/bin/bash
#SBATCH --job-name=build_squakr_initial_index     # Job name
#SBATCH --output=./logs/%x-%j.out    # Standard output file
#SBATCH --error=./logs/%x-%j.err    # Standard error file
#SBATCH --partition=short     # Partition/queue name
#SBATCH --nodes=1              # Number of nodes/machines
#SBATCH --ntasks=1             # Number of tasks/separate processes
#SBATCH --cpus-per-task=16      # CPU cores per task
#SBATCH --mem=32G               # amount of ram
#SBATCH --time=05:00:00        # Time limit hrs:min:sec

## load modules ##
module purge
module load Boost/1.88.0
## run vars ##
export kmer_size=28
export slots=31
export write_dir="/scratch/w.galbraith/CS7800_group_4/mantis/mantis_output"
export fastq_dir="/scratch/w.galbraith/CS7800_group_4/mantis/sra_data/toy_example_files" ## dir with input fastqs ##
export threads=16 ## make sure this matches the sbatch config

## bins/lib paths ##
export SQUEAKR_BIN="/scratch/w.galbraith/CS7800_group_4/mantis/bin/squeakr"
export LD_LIBRARY_PATH="/scratch/w.galbraith/CS7800_group_4/mantis/lib/boost-1.84.0/"
## run method ##
export mb_size=100000 ## 100 mb size
export squeakr_dir="${write_dir}/squeakr_files_mini"
export index_dir="${write_dir}/index"
mkdir -p $squeakr_dir $index_dir
## get squeakr files ##
echo "Creating squeakr files..."
for fastq_file in $fastq_dir/*.fastq; do
	## determine cutoff ## 
        size=$(du -k "$fastq_file" | awk -F '\t' '{print $1 }')
        if ((3* mb_size >= size)); then
                export cutoff=1;
        elif ((5 * mb_size >= size)); then
                export cutoff=3;
        elif ((10 * mb_size >= size)); then
                export cutoff=10;
        elif ((30 * mb_size >= size)); then
                export cutoff=20;
        else
                export cutoff=50;
        fi
        base=$(basename "$fastq_file" .fastq)
        write_name=$squeakr_dir/$base.squeakr
        echo $write_name
        $SQUEAKR_BIN count \
		-e \
                --no-counts \
                --cutoff $cutoff \
                -k $kmer_size \
                -s $slots \
                -t $threads \
                -o $write_name \
                $fastq_file
done
echo "squeakr done"
