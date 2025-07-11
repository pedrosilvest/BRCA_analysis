#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project
script_name=FilterBam
script_number=06
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



# Define the path to the Singularity container image with Samtools
SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"


file_paths=($(printf '"%s"\n' "${file_paths[@]}"))

INPUT="\${file_paths[\$SLURM_ARRAY_TASK_ID]}"
echo "\$INPUT"

sample=\$(basename "\$INPUT")

INPUT=(\$INPUT/bam/\${sample}.bam)

echo "Input file: \$INPUT"
echo "Dest dir: \$(dirname "\$INPUT")"
# Extract the sample name from the input BAM file
sample=\$(basename "\$INPUT" .bam)
echo "Sample: \$sample"

# Create a temporary output directory
temp_output="$temp_path/\${sample}_filter_\${SLURM_JOB_ID}"
mkdir -p "\$temp_output"

# Print debugging information
echo "Processing BAM file: \$INPUT"
echo "Output directory: \$temp_output"

# Filtering
echo "FILTERING"
filtered_bam="\$temp_output/\$(basename "\$INPUT")"

srun apptainer exec "\$SAMTOOLS_IMAGE" samtools view -@ 30 -q 10 -F 256 -F 4 -b "\$INPUT" > "\$filtered_bam"

# Quick check on the filtered BAM file
if srun apptainer exec "\$SAMTOOLS_IMAGE" samtools quickcheck "\$filtered_bam"; then
    echo "The BAM file \$filtered_bam passed the quick check and is not corrupted."
    
    # Sorting
    echo "SORTING"
    sorted_file="\$temp_output/\${sample}_sorted.bam"
    srun apptainer exec "\$SAMTOOLS_IMAGE" samtools sort -@ 30 -m 4G "\$filtered_bam" -o "\$sorted_file"
    
    # Indexing
    echo "INDEXING"
    srun apptainer exec "\$SAMTOOLS_IMAGE" samtools index -@ 30 "\$sorted_file"
    
    # Stats
    echo "GENERATING STATS"
    STATS_OUTPUT="\$temp_output/\${sample}_sorted_stats.txt"
    srun apptainer exec "\$SAMTOOLS_IMAGE" samtools flagstat -@ 30 "\$sorted_file" > "\$STATS_OUTPUT"

fi


BAM_DIR=$output_path/\$sample/bam
STATS_DIR=$output_path/\$sample/stats

echo "MOVING \$sorted_file to \$BAM_DIR"
mv "\$sorted_file" "\$BAM_DIR"
mv "\$sorted_file.bai" "\$BAM_DIR"
mv "\$STATS_OUTPUT" "\$STATS_DIR"

rm -rf \$temp_output





EOL


