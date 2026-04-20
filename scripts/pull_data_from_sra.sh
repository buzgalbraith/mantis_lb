#!/bin/bash
#SBATCH --job-name=pull_sra_data     # Job name
#SBATCH --output=./logs/%x-%j.out    # Standard output file
#SBATCH --error=./logs/%x-%j.err    # Standard error file
#SBATCH --partition=short     # Partition/queue name
#SBATCH --nodes=1              # Number of nodes/machines
#SBATCH --ntasks=1             # Number of tasks/separate processes
#SBATCH --cpus-per-task=16      # CPU cores per task
#SBATCH --mem=64G               # amount of ram
#SBATCH --time=12:00:00        # Time limit hrs:min:sec
## load modules ## 
module load sratoolkit/12Dec2024

## run args ## 
export sra_dir="/scratch/w.galbraith/CS7800_group_4/mantis/sra_data/"
export files_to_pull=750
export initial_index_files=500
export load_balance_files=250
export n_threads=16
## run method ## 
export fastq_files_dir="${sra_dir}/fastq_files"
export prefetch_files_dir="${sra_dir}/prefetch_files"
export load_ballance_files_dir="${sra_dir}/load_ballance_fastq_files"
export initial_index_files_dir="${sra_dir}/initial_index_fastq_files" 
mkdir -p $fastq_files_dir  $initial_index_files_dir $load_ballance_files_dir $prefetch_files_dir
curl https://www.cs.cmu.edu/~ckingsf/software/bloomtree/srr-list.txt | awk '{print $2}' > ${sra_dir}/accessions.txt

tail -n "$files_to_pull" ${sra_dir}/accessions.txt | while read acc; do
    echo "Fetching $acc..."
    prefetch "$acc" --output-directory $prefetch_files_dir
    fasterq-dump $prefetch_files_dir/"$acc"/"$acc".sra \
        --outdir $fastq_files_dir \
        --threads $n_threads \
        --split-files
done

echo "Saving initial index files to ${initial_index_files_dir}"
ls $fastq_files_dir | head -n $initial_index_files | xargs -I {} mv $fastq_files_dir/{} $initial_index_files_dir/{}

echo "Saving load ballance files to ${load_ballance_files_dir}"
ls $fastq_files_dir | head -n $load_balance_files | xargs -I {} mv $fastq_files_dir/{} $load_ballance_files_dir/{}

