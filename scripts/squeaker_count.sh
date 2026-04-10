#!/bin/bash
#SBATCH --job-name=build_squakr      # Job name
#SBATCH --output=./logs/%x-%j.out    # Standard output file
#SBATCH --error=./logs/%x-%j.err    # Standard error file
#SBATCH --partition=sharing     # Partition/queue name
#SBATCH --nodes=1              # Number of nodes/machines
#SBATCH --ntasks=1             # Number of tasks/separate processes
#SBATCH --cpus-per-task=16      # CPU cores per task
#SBATCH --mem=32G               # amount of ram
#SBATCH --time=01:00:00        # Time limit hrs:min:sec

## load modules ##
module purge
module load Boost/1.88.0
## run vars ##
export kmer_size=28
export slots=31
export write_dir="./mantis_output"
export fastq_dir="/scratch/w.galbraith/CS7800_group_4/mantis/sra_data/toy_example_files" ## dir with input fastqs ##
export threads=16 ## make sure this matches the sbatch config

## bins/lib paths ##
export SQUEAKR_BIN="/scratch/w.galbraith/CS7800_group_4/mantis/bin/squeakr"
export LD_LIBRARY_PATH="/scratch/w.galbraith/CS7800_group_4/mantis/lib/boost-1.84.0/"
## run method ##
export squeakr_dir="${write_dir}/squeakr_files"
export index_dir="${write_dir}/index"
mkdir -p $squeakr_dir $index_dir
## get squeakr files ##
echo "Creating squeakr files..."
for fastq_file in $fastq_dir/*.fastq; do
        base=$(basename "$fastq_file" .fastq)
        write_name=$squeakr_dir/$base.squeakr
        echo $write_name
        $SQUEAKR_BIN count \
                --no-counts \
                --cutoff 3 \
                -k $kmer_size \
                -s $slots \
                -t $threads \
                -o $write_name \
                $fastq_file
done
echo "squeakr done"
