exp_mode=${1:-1}
mean_mode=${2:-1}
exp_dir='/scratch/w.galbraith/expirements/'
header=1
if [[ $exp_mode == "1" ]]; then
	header=1
	for i in 250 500 1000 1500; do
		res_file="${exp_dir}/${i}_sample_ids_and_5_clusters/results/summary/mean_variance_change_summary.txt"
		python expirements/utils/markdown_table.py $res_file $header --experiments $mean_mode
		header=0
	done;
else
for i in 3 5 7 9; do
		res_file="${exp_dir}/1500_sample_ids_and_${i}_clusters/results/summary/mean_variance_change_summary.txt"
		python expirements/utils/markdown_table.py $res_file $header --workers $mean_mode
		header=0
	done;
fi
