#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project
script_name=MergeBams
script_number=04
mem=240G
cpus_per_task=40


#create folders
logs_path=$scratch_path/output/logs
scripts_path=$scratch_path/scripts
temp_path=$scratch_path/temp

#mkdir -p $output_path
mkdir -p $logs_path
mkdir -p $scripts_path
mkdir -p $temp_path

# Should be a list of paths!!
file_paths=(
"/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/CBE7negB"
"/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/CBE7ttdB"
)

# Get the total number of files (array length)
total_files=${#file_paths[@]}

# Subtract 1 from the total number of files
total_arrays=$((total_files - 1))


setup_sbatch="$scripts_path/${script_number}_${script_name}.sbatch"

echo "Creating scripts... $setup_sbatch"
cat >"$setup_sbatch" <<EOL
#!/bin/bash
#SBATCH --job-name=$script_name
#SBATCH --output=$logs_path/${script_number}_${script_name}_%A_%a.out
#SBATCH --array=0-$total_arrays%5
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --mem=$mem
#SBATCH --cpus-per-task=$cpus_per_task


# MergeBams and Get Coverage stats

file_paths=($(printf '"%s"\n' "${file_paths[@]}"))

INPUT="\${file_paths[\$SLURM_ARRAY_TASK_ID]}"
echo "\$INPUT"

bam_files=(\$INPUT/bam/*_sorted.bam)

#echo \${bam_files[@]}

sample=\$(basename "\$INPUT")


temp_output="/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/temp/\${sample}_\${SLURM_JOB_ID}"
mkdir -p "\$temp_output"


SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"


output_bam="\$temp_output/\${sample}.bam"
echo "merging bams"
srun apptainer exec \$SAMTOOLS_IMAGE samtools merge -@ 30 \$output_bam "\${bam_files[@]}"

echo "indexing"
srun apptainer exec \$SAMTOOLS_IMAGE samtools index -@ 30 \$output_bam

# Run

output_coverage="\$temp_output/\${sample}_coverage.txt"

echo "samtools depth to get coverage"
srun apptainer exec "\$SAMTOOLS_IMAGE" samtools depth \$output_bam > \$output_coverage


#-a to include 0 coverage positions
echo "calculating average coverage"
awk '{sum+=\$3} END {print "Average coverage: " sum/NR}' \$output_coverage


# STATISTICS
stats_output="\$temp_output/\${sample}_stats.txt"

echo "samtools flagstats"
srun apptainer exec "\$SAMTOOLS_IMAGE" samtools flagstat -@ 30 \$output_bam > \$stats_output

echo "RG:"
srun apptainer exec "\$SAMTOOLS_IMAGE" samtools view -H \$output_bam | grep "^@RG"


echo "moving files"

mv \$output_bam* \$INPUT/bam
mv \$output_coverage \$INPUT/stats
mv \$stats_output \$INPUT/stats

rm -rf \$temp_output


EOL