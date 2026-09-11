# Mantis_lb (Mantis load balancer)
Mantis_lb is a load balancing framework for [Mantis](https://github.com/splatlab/mantis) 
sequence search indices. It partitions FASTQ files into sub-indices using either 
round-robin or sketch-based similarity clustering, enabling scalable colored de Bruijn 
graph construction across large SRA datasets
## Package structure 
- `src/` Rust code for load balancing
- `scripts/` Code for experiment evaluation and building Mantis sub-indices

## Dependencies
- Mantis: version `0.2.0` (commit [0fb7dbb](https://github.com/splatlab/mantis/tree/0fb7dbb60e4a38aa21da7c85f98565786af5fe62))
- Squeakr: version `0.7` (commit [dcfaa18](https://github.com/splatlab/squeakr/tree/dcfaa18f267814d9e7d3437fbfc7348b869dab88)) 

## Reproduce results
To reproduce the results run the following scripts. Be sure to change the run arguments at the top of each script as required.
## Pulling data and set up
1. Pull the data from SRA tools and split it into initial index and load balance files with `sbatch scripts/pull_data_from_sra.sh`
1. Use Squeakr to build CQFs for all of FASTQ files `scripts/build_cqfs.sh`
## Running experiments
The bellow will need to be run for each experiment. Axis to chose are how many worker nodes to have, have many files to distribute and how many files to use in the initial index
1. Split the FASTQ files randomly by sample ID into initial index and files to distribute. Run `bash scripts/split_sample_ids.sh`
1. Build sketches of the initial worker content `sbatch scripts/build_initial_worker_sketches.sh`
1. Distribute the files using each load balancing method `sbatch scripts/run_load_ballance.sh`
1. Merge the initial worker index assignment with the load balancing index assignment `bash scripts/concat_index_assigment.sh`
1. Build the Mantis indexes `sbatch scripts/build_mantis_index.sh`. **note after this runs check the log to make sure all of the indexes built correctly**
## Getting results ## 
The bellow are some utility files i built for quickly gathering results
1. Can gather raw results with `bash expirements/gather_results.sh <jid>` where `<jid>` is the job id of your mantis run
1. Can generate a markdown table of mean and variance summary. `bash scripts/generate_markdown_results.sh> <exp> <mean>`  where `<exp>` indicates if the table should be across experiments (otherwise across workers), and `<mean>` denotes if we are looking at change in mean or variance
