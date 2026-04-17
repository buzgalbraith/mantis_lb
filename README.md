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
1. Pull the data from SRA tools and split it into initial index and load balance files with `sbatch scripts/pull_data_from_sra.sh`
1. Build sketches of the initial index `sbatch scripts/build_initial_sketch_clusters.sh`
1. Do round-robin based cluster assignment `sbatch scripts/run_round_robin_cluster.sh` (this should only take a second)
1. Do similarity based cluster assignment `sbatch scripts/run_similarity_cluster.sh`
1. Use Squeakr to get CQFs for all files you want to load balance `sbatch scripts/squeaker_count_load_ballance.sh`
1. Create the Mantis sub-indices with Round-robin assignment `sbatch scripts/build_mantis_index_round_robin.sh`
1. Create the Mantis sub-indices with similarity assignment `sbatch scripts/build_mantis_index_similarity.sh`