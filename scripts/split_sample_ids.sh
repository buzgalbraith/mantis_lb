export fastq_dir='/projects/gyorilab/buz/mantis_lb/fastqs'
export fastq_list_dir='/scratch/w.galbraith/expirements/fastq_list/' 
export initial_index_number=250
export load_ballance_number=750

mkdir -p $fastq_list_dir
## get unique sample ids that were pulled (ie group paired reads) then shuffle them ##
find $fastq_dir -name "*.fastq.gz" | xargs -I {} basename {} | \
	sed 's/\.fastq.gz//' | sed 's/_[12]//' | \
	sort -u | shuf \
	> "${fastq_dir}/sample_ids.txt"
## now split those and get a list of associated files for each sample id
head -n $initial_index_number "${fastq_dir}/sample_ids.txt" | xargs -I {} sh -c 'ls '"$fastq_dir"'/{}*'> "${fastq_list_dir}/initial_index_files_from_${initial_index_number}_sample_ids.txt"
tail -n +$((initial_index_number+1)) "${fastq_dir}/sample_ids.txt" | head -n $load_ballance_number | xargs -I {} sh -c 'ls '"$fastq_dir"'/{}*' > "${fastq_list_dir}/load_balance_files_from_${load_ballance_number}_sample_ids.txt"
