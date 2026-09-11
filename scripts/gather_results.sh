#!/bin/bash

job_id=$1
export base_expirements="/scratch/w.galbraith/expirements"
export n_clusters=5
export num_lb_sample_ids=500 # number of sample ids for total index
expirement_dir="${base_expirements}/${num_lb_sample_ids}_sample_ids_and_${n_clusters}_clusters"
results_dir="${expirement_dir}/results"
kmer_dir="${results_dir}/kmer_counts"
color_class_dir="${results_dir}/color_class_counts"
file_count_dir="${results_dir}/file_counts"
size_dir="${results_dir}/index_size"
summary_dir="${results_dir}/summary"
rm -rf $results_dir # remove existing summary
mkdir -p $kmer_dir $color_class_dir $size_dir $summary_dir $file_count_dir
declare -A index_to_lb_method=(
    [0]="round_robin"
    [1]="weighted_random"
    [2]="similarity"
)
idx=0
kmer_summary="${summary_dir}/kmer_summary.txt"
echo "method kmer_count mean variance" > $kmer_summary
color_class_summary="${summary_dir}/color_class_summary.txt"
echo "method color_class_count mean variance" > $color_class_summary
file_counts_summary="${summary_dir}/file_counts_summary.txt"
echo "method file_count mean variance" > $file_counts_summary
overall_index_size_summary="${summary_dir}/overall_index_size_summary.txt"
echo "method overall_index_size mean variance" > $overall_index_size_summary
cqf_size_summary="${summary_dir}/cqf_size_summary.txt"
echo "method cqf_size mean variance" > $cqf_size_summary
ct_size_summary="${summary_dir}/color_table_size_summary.txt"
echo "method ct_size mean variance" > $ct_size_summary
mst_ct_index_size_summary="${summary_dir}/mst_color_table_size_summary.txt"
echo "method mst_ct_size mean variance" > $mst_ct_index_size_summary
run_time_summary="${summary_dir}/run_time.txt"
echo "method run_time" > $run_time_summary
mean_variance_change_summary="${summary_dir}/mean_variance_change_summary.txt"
echo "metric mean_change variance_change" > $mean_variance_change_summary
echo "Writing results to ${results_dir}"
echo "Summary available in ${summary_dir}"
for log in ./logs/*${job_id}*.out; do
    echo $log
    method="${index_to_lb_method[${idx}]}"
    kmer_path="${kmer_dir}/${method}_kmer_counts.txt"
    cat $log | grep 'Final colored dBG has' | awk '{ print $9 }' > $kmer_path
    mv=$(python ./expirements/utils/get_stats.py $kmer_path)
    cat $log | grep 'Final colored dBG has' | awk '{ print $9 }' | paste -sd+ | bc | xargs -I {} echo $method {} $mv >> $kmer_summary
    color_class_path="${color_class_dir}/${method}_color_class_counts.txt"
    cat $log | grep 'Final colored dBG has' | awk '{ print $12 }' > $color_class_path
    mv=$(python ./expirements/utils/get_stats.py $color_class_path)
    cat $log | grep 'Final colored dBG has' | awk '{ print $12 }' | paste -sd+ | bc | xargs -I {} echo $method {} $mv >> $color_class_summary
    file_count_path="${file_count_dir}/${method}_file_counts.txt"
    cat $log | grep '# of experiments:' | awk '{ print $8 }' >> $file_count_path
    mv=$(python ./expirements/utils/get_stats.py $file_count_path)
    cat $log | grep '# of experiments:' | awk '{ print $8 }' | paste -sd+ | bc | xargs -I {} echo $method {} $mv >> $file_counts_summary
    overall_size_path="${size_dir}/${method}_overall_index_size.txt"
    cqf_size_path="${size_dir}/${method}_cqf_size.txt"
    ct_size_path="${size_dir}/${method}_color_table_size.txt"
    mst_ct_size_path="${size_dir}/${method}_mst_color_table_size.txt"
    for i in $(seq 0 $(($n_clusters-1))); do
        index_dir="${expirement_dir}/${method}_index/cluster_${i}_index"
        du $index_dir/*.ser | awk '{print $1 }'| paste -sd+ | bc | xargs -I {} echo {} >> $cqf_size_path ## MST CT
        du $index_dir/*.cls | awk '{print $1 }'| paste -sd+ | bc | xargs -I {} echo {} >> $ct_size_path ## MST CT
        du $index_dir/*.bv | awk '{print $1 }'| paste -sd+ | bc | xargs -I {} echo {} >> $mst_ct_size_path ## MST CT
        size=$(du $index_dir | awk '{print $1 }')
        echo $size >> $overall_size_path        
    done
    mv=$(python ./expirements/utils/get_stats.py $overall_size_path)
    awk '{ sum += $1 } END { print sum }' $overall_size_path | xargs -I {} echo $method {} $mv >> $overall_index_size_summary
    mv=$(python ./expirements/utils/get_stats.py $cqf_size_path)
    awk '{ sum += $1 } END { print sum }' $cqf_size_path | xargs -I {} echo $method {} $mv >> $cqf_size_summary
    mv=$(python ./expirements/utils/get_stats.py $ct_size_path)
    awk '{ sum += $1 } END { print sum }' $ct_size_path | xargs -I {} echo $method {} $mv >> $ct_size_summary
    mv=$(python ./expirements/utils/get_stats.py $mst_ct_size_path)
    awk '{ sum += $1 } END { print sum }' $mst_ct_size_path | xargs -I {} echo $method {} $mv >> $mst_ct_index_size_summary
    seff "${job_id}_${idx}" | grep "Job Wall-clock time" | awk ' { print $4 }' | xargs -I {} echo $method {} >> $run_time_summary
    idx=$((idx+1))
done
python ./expirements/utils/get_change.py $kmer_summary | xargs -I {} echo kmer_count {} >> $mean_variance_change_summary
python ./expirements/utils/get_change.py $color_class_summary | xargs -I {} echo color_class_count {} >> $mean_variance_change_summary
python ./expirements/utils/get_change.py $file_counts_summary | xargs -I {} echo file_count {} >> $mean_variance_change_summary
python ./expirements/utils/get_change.py $overall_index_size_summary | xargs -I {} echo overall_index_size {} >> $mean_variance_change_summary
python ./expirements/utils/get_change.py $cqf_size_summary | xargs -I {} echo cqf_size {} >> $mean_variance_change_summary
python ./expirements/utils/get_change.py $ct_size_summary | xargs -I {} echo color_class_table_size {} >> $mean_variance_change_summary
python ./expirements/utils/get_change.py $mst_ct_index_size_summary | xargs -I {} echo mst_color_class_table_size {} >> $mean_variance_change_summary
