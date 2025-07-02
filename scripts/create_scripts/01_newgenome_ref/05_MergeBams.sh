#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples/newref
script_name=MergeBams
script_number=05
mem=240G
cpus_per_task=40

# IMAGES
BWA="/mnt/beegfs/apptainer/images/bwa_latest.sif"
SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"

GENOME="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples/reference/CBE7negB_reference.fa"

#create folders
logs_path=$scratch_path/output/logs/01_newgenome_ref
scripts_path=$scratch_path/scripts/01_newgenome_ref
temp_path=$scratch_path/temp

#mkdir -p $output_path
mkdir -p $logs_path
mkdir -p $scripts_path
mkdir -p $temp_path


total_arrays=$(find $output_path -mindepth 1 -maxdepth 1 -type d | wc -l)

total_arrays=$((total_arrays - 1))



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

file_paths=($(find $output_path -mindepth 1 -maxdepth 1 -type d))

INPUT="\${file_paths[\$SLURM_ARRAY_TASK_ID]}"

echo "\$INPUT"

bam_files=(\$INPUT/bam/*_sorted.bam)

#echo \${bam_files[@]}

sample=\$(basename "\$INPUT")


temp_output="$temp_path/\${sample}_\${SLURM_JOB_ID}"
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

mv \$output_bam* $output_path/\$sample/bam
mv \$output_coverage $output_path/\$sample/stats
mv \$stats_output $output_path/\$sample/stats

rm -rf \$temp_output


EOL