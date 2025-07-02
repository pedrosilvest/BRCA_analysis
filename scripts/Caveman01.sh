normal_bam=(/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/C2040neg/bam/C2040neg_sorted.dup.bqsr.bam
            /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/C2041neg/bam/C2041neg_sorted.dup.bqsr.bam)

tumor_bam=(/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/C2040ttd/bam/C2040ttd_sorted.dup.bqsr.bam
            /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/C2041ttd/bam/C2041ttd_sorted.dup.bqsr.bam)

sbatch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/scripts/caveman
mkdir -p $sbatch_path

for i in "${!normal_bam[@]}"; do


    reference_index="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/unzipped/GRCh38.primary.genome.fa.fai"
    ignore_regions="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/ignore_regions.tsv"
    >"$ignore_regions"

    temp_path="/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/temp/C204${i}_caveman"
    log_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/output/logs/caveman/C204${i}
    mkdir -p $log_path

    other_path="/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/temp/C204${i}_caveman/C204${i}_caveman_others"
    results_path="/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/temp/C204${i}_caveman/C204${i}_caveman_results"
    error_path="/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/temp/C204${i}_caveman/C204${i}_caveman_error"
    mkdir -p "$temp_path"
    mkdir -p "$other_path"
    mkdir -p "$results_path"
    mkdir -p "$error_path"

    config_file="$other_path/C204${i}_caveman.cfg.ini"
    split_file="$other_path/splitList"
    alg_bean_file="$other_path/alg_bean"
    covs_arr="/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/temp/C204${i}_caveman/C204${i}_caveman_others/C204${i}_caveman.covs_arr"
    probs_arr="/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/temp/C204${i}_caveman/C204${i}_caveman_others/C204${i}_caveman.probs_arr"
  

    caveman_image=/mnt/beegfs/apptainer/images/caveman_latest.sif

    echo "Creating scripts..."
    setup_sbatch="$sbatch_path/01_Caveman_C204${i}_setup.sbatch"

    cat >"$setup_sbatch" <<EOL
#!/bin/bash
#SBATCH --job-name=cav_setup
#SBATCH --output=$log_path/caveman01_setup_%A_%a.out
#SBATCH --time=72:00:00
#SBATCH --array=0
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=40

caveman=$caveman_image

normal_bam="${normal_bam[i]}"
tumoral_bam="${tumor_bam[i]}"

apptainer exec "\$caveman" caveman setup \\
    -n "\$normal_bam" \\
    -t "\$tumoral_bam" \\
    -r "$reference_index" \\
    -g "$ignore_regions" \\
    -f "$results_path" \\
    -c "$config_file" \\
    -a "$alg_bean_file" \\
    -l "$split_file"

echo "setup created"

EOL


echo "split caveman sbatch file"

    total_arrays=$(cat "$reference_index" | wc -l)
    setup_sbatch="$sbatch_path/02_Caveman_C204${i}_split.sbatch"

    cat >"$setup_sbatch" <<EOL
#!/bin/bash
#SBATCH --job-name=cav_split
#SBATCH --output=$log_path/caveman02_split_%A_%a.out
#SBATCH --time=72:00:00
#SBATCH --array=1-$total_arrays%5
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=40

caveman=$caveman_image

index=\$SLURM_ARRAY_TASK_ID

echo "doing \$index for ${normal_bam[i]}"
apptainer exec \$caveman caveman split \\
    -i \$index \\
    -f $config_file

EOL

    setup_sbatch="$sbatch_path/03_Caveman_C204${i}_mstep.sbatch"

    cat >"$setup_sbatch" <<EOL
#!/bin/bash
#SBATCH --job-name=cav_mstep
#SBATCH --output=$log_path/caveman03_mstep_%A_%a.out
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=40
#SBATCH --array=0-29%4


sample=C204${i}
files=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/temp/\${sample}_caveman/\${sample}_caveman_others/splitList*
split_file=$split_file

cat \$files > \$split_file

caveman=$caveman_image

total_arrays=\$(cat "$split_file" | wc -l)


# Calculate the size of each interval

number_intervals=30

interval_size=\$((\$total_arrays / \$number_intervals))

# Calculate the remainder to distribute among the intervals
remainder=\$((\$total_arrays  % \$number_intervals))

# Initialize the starting point of the intervals
start=0

# Initialize an empty array to store the intervals
intervals=()

# Loop through each interval
for ((i = 0; i < \$number_intervals; i++)); do
    # Calculate the end point of the interval
    end=\$((start + interval_size))
    
    # If there's a remainder, distribute it among the intervals
    if [ \$remainder -gt 0 ]; then
        ((end++))
        ((remainder--))
    fi
    
    # Add the interval to the array
    intervals+=("\$start-\$end")
    
    # Update the starting point for the next interval
    start=\$end
done

# Print the intervals stored in the array
echo "\${intervals[@]}"


interval=\${intervals[\$SLURM_ARRAY_TASK_ID]}
start_interval=\$(echo "\$interval" | cut -d'-' -f1)
end_interval=\$(echo "\$interval" | cut -d'-' -f2)

error_file=$error_path/mstep_error_indexes.txt

# Perform the caveman for each index in the interval
for ((index = \$start_interval; index <= \$end_interval; index++)); do

    srun apptainer exec "\$caveman" caveman mstep \
        -i "\$index" \
        -f "$config_file"

    # Check if there was an error
    if [ $? -ne 0 ]; then
        echo "Error occurred for index \$index: \$output"
        # Save the index to the error file
        echo \$index >> \$error_file
    fi

done


EOL
 
   setup_sbatch="$sbatch_path/04_Caveman_C204${i}_merge.sbatch"

    cat >"$setup_sbatch" <<EOL
#!/bin/bash
#SBATCH --job-name=cav_merge
#SBATCH --output=$log_path/caveman04_merge_%A_%a.out
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --array=1

caveman=$caveman_image

apptainer exec "\$caveman" caveman merge \\
    -f "$config_file" \\
    -c "$covs_arr" \\
    -p "$probs_arr"


EOL

   setup_sbatch="$sbatch_path/05_Caveman_C204${i}_estep.sbatch"

    cat >"$setup_sbatch" <<EOL
#!/bin/bash
#SBATCH --job-name=cav_estep
#SBATCH --output=$log_path/caveman05_estep_%A_%a.out
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=40
#SBATCH --array=0-29%8

caveman=$caveman_image

total_arrays=\$(cat "$split_file" | wc -l)


# Calculate the size of each interval

number_intervals=30

interval_size=\$((\$total_arrays / \$number_intervals))

# Calculate the remainder to distribute among the intervals
remainder=\$((\$total_arrays  % \$number_intervals))

# Initialize the starting point of the intervals
start=0

# Initialize an empty array to store the intervals
intervals=()

# Loop through each interval
for ((i = 0; i < \$number_intervals; i++)); do
    # Calculate the end point of the interval
    end=\$((start + interval_size))
    
    # If there's a remainder, distribute it among the intervals
    if [ \$remainder -gt 0 ]; then
        ((end++))
        ((remainder--))
    fi
    
    # Add the interval to the array
    intervals+=("\$start-\$end")
    
    # Update the starting point for the next interval
    start=\$end
done

# Print the intervals stored in the array
echo "\${intervals[@]}"



interval=\${intervals[\$SLURM_ARRAY_TASK_ID]}
start_interval=\$(echo "\$interval" | cut -d'-' -f1)
end_interval=\$(echo "\$interval" | cut -d'-' -f2)
    
# Perform the caveman for each index in the interval
for ((index = \$start_interval + 1; index <= \$end_interval; index++)); do

    srun apptainer exec "\$caveman" caveman estep \\
        -i "\$index" \\
        -f "$config_file" \\
        -g "$covs_arr" \\
        -o "$probs_arr" \\
        -v GRCh38 \\
        -w Human

done



EOL


setup_sbatch="$sbatch_path/06_Gatk_C204${i}_merge.sbatch"

    cat >"$setup_sbatch" <<EOL
#!/bin/bash
#SBATCH --job-name=gatk_merge
#SBATCH --output=$log_path/gatk_merge_%A_%a.out
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40
#SBATCH --array=0


######### IMAGE #########
GATK_IMAGE="/mnt/beegfs/apptainer/images/gatk4.sif"
SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"
########


input_directory="$results_path"

final_dir="$results_path/0_all_vcf"

temp_output="/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/temp/mergeVCF_caveman_\${SLURM_JOB_ID}"

# Create the output directory if it doesn't exist
mkdir -p "\$temp_output"

find \$input_directory -name "*.muts.vcf.gz" > \$temp_output/input_vcf_files.list

srun apptainer exec "\$GATK_IMAGE" gatk MergeVcfs \\
    -I \$temp_output/input_vcf_files.list \\
    -O \$temp_output/all_muts.vcf.gz


find \$input_directory -name "*.snps.vcf.gz" > \$temp_output/input_vcf_files.list

srun apptainer exec "\$GATK_IMAGE" gatk MergeVcfs \\
    -I \$temp_output/input_vcf_files.list \\
    -O \$temp_output/all_snps.vcf.gz


mkdir -p \$final_dir

mv \$temp_output/* \$final_dir

copy $results_path/*/*.no_analysis.bed \$final_dir/no_analysis.bed

rm -rf \$temp_output

EOL

done