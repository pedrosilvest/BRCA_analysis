#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project
script_name=MarkDuplicates
script_number=07
mem=240G
cpus_per_task=40



GENOME="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/GRCh38.primary.genome.fa.gz"

#create folders
logs_path=$scratch_path/output/logs
scripts_path=$scratch_path/scripts
temp_path=$scratch_path/temp

#mkdir -p $output_path
mkdir -p $logs_path
mkdir -p $scripts_path
mkdir -p $temp_path


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

file_paths=($(printf '"%s"\n' "${file_paths[@]}"))

INPUT="\${file_paths[\$SLURM_ARRAY_TASK_ID]}"
echo "\$INPUT"

sample=\$(basename "\$INPUT")

INPUT=(\$INPUT/bam/\${sample}_sorted.bam)

input_file=\$INPUT
echo "Input file: \$input_file"
echo "Sample: \$sample"



# Define the output directory
temp_output="$temp_path/\${sample}_markdup_\${SLURM_JOB_ID}"

# Create the output directory if it doesn't exist
mkdir -p "\$temp_output"


################## output_file name is .dup.bam ############### 

output_file=\$temp_output/\$(basename "\$input_file" .bam).dup.bam
metrics_file=\$temp_output/\$(basename "\$input_file" .bam)_metrics_dup.txt

echo "Output: \$output_file"
echo "Metrics: \$metrics_file"

######### IMAGE #########
GATK_IMAGE="/mnt/beegfs/apptainer/images/gatk_latest.sif"
SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"
########


srun apptainer exec "\$SAMTOOLS_IMAGE" samtools quickcheck "\$input_file"

if [ \$? -eq 0 ]; then
    echo "The BAM file \$input_file is ok"
else
    echo "The BAM file \$input_file failed the quick check and may be corrupted."
fi


START_TIME=\$(date +%s)

srun apptainer exec "\$GATK_IMAGE" gatk --java-options "-Xmx120G -XX:ParallelGCThreads=36 -XX:ConcGCThreads=36" MarkDuplicates \
  -I "\$input_file" \
  -O "\$output_file" \
  -M "\$metrics_file" \
  --TMP_DIR "\$temp_output"


END_TIME=\$(date +%s)
ELAPSED_TIME=\$((END_TIME - START_TIME))

echo "Job finished in \$((\$ELAPSED_TIME / 3600)) hours, \$(((\$ELAPSED_TIME % 3600) / 60)) minutes, and \$((\$ELAPSED_TIME % 60)) seconds."


output_file_DIR=\$(dirname "\$INPUT")

metrics_file_DIR=\$(dirname "\$INPUT" | sed 's/bam\$/markdup_metrics/')
mkdir -p \$metrics_file_DIR


# Indexing
echo "INDEXING"
srun apptainer exec "\$SAMTOOLS_IMAGE" samtools index -@ 30 "\$output_file"
    
if [ \$? -eq 0 ]; then
    mv \$output_file* \$output_file_DIR
    mv \$metrics_file \$metrics_file_DIR
    rm -rf \$temp_output
else
    echo "The BAM file \$input_file failed the quick check and may be corrupted."
fi


EOL