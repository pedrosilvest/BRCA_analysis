#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples/newref
script_name=FilterBam
script_number=06
mem=240G
cpus_per_task=40

# IMAGE

SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"


#create folders
logs_path=$scratch_path/output/logs/01_newgenome_ref
scripts_path=$scratch_path/scripts/01_newgenome_ref
temp_path=$scratch_path/temp

#mkdir -p $output_path
mkdir -p $logs_path
mkdir -p $scripts_path
mkdir -p $temp_path


total_arrays=$(find $output_path/*/bam -type f -name "*.bam" | wc -l)
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

####################
SAMTOOLS_IMAGE=$SAMTOOLS_IMAGE
####################

# Define an array of file paths for BAM files to process
file_path=\$(find $output_path/*/bam -type f -name "*.bam")


echo "all paths to do: \$file_path"

IFS=$'\r\n' GLOBIGNORE='*' command eval 'INPUTS=(\$(echo "\$file_path"))'

INPUT="\${INPUTS[\$SLURM_ARRAY_TASK_ID]}"

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
    srun apptainer exec "\$SAMTOOLS_IMAGE" samtools sort -@ 30 "\$filtered_bam" -o "\$sorted_file"
    
    # Indexing
    echo "INDEXING"
    srun apptainer exec "\$SAMTOOLS_IMAGE" samtools index -@ 30 "\$sorted_file"
    
    # Stats
    echo "GENERATING STATS"
    STATS_OUTPUT="\$temp_output/\${sample}_sorted_stats.txt"
    srun apptainer exec "\$SAMTOOLS_IMAGE" samtools flagstat -@ 30 "\$sorted_file" > "\$STATS_OUTPUT"

fi


BAM_DIR=\$(dirname "\$INPUT")
STATS_DIR=$output_path/\$sample/stats

echo "MOVING \$sorted_file to \$BAM_DIR"
mv "\$sorted_file" "\$BAM_DIR"
mv "\$sorted_file.bai" "\$BAM_DIR"
mv "\$STATS_OUTPUT" "\$STATS_DIR"

rm -rf \$temp_output





EOL