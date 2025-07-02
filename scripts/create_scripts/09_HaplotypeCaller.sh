#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples
script_name=HaplotypeCaller
script_number=09
mem=240G
cpus_per_task=40

# arrays setup
total_arrays=$(find $output_path/*/bam -type f -name "*bqsr.bam" | wc -l)
total_arrays=$((total_arrays - 1))


GENOME="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/GRCh38.primary.genome.fa.gz"

#create folders
logs_path=$scratch_path/output/logs
scripts_path=$scratch_path/scripts
temp_path=$scratch_path/temp

#mkdir -p $output_path
mkdir -p $logs_path
mkdir -p $scripts_path
mkdir -p $temp_path





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

file_path=\$(find $output_path/*/bam -type f -name "*_sorted.dup.bqsr.bam")

echo "all file paths: \$file_path"

IFS=\$'\r\n' GLOBIGNORE='*' command eval 'INPUTS=(\$(echo "\$file_path"))'


INPUT="\${INPUTS[\$SLURM_ARRAY_TASK_ID]}"

sample=\$(basename "\$INPUT" _sorted.dup.bqsr.bam)



# Construct the full path to the input BAM file
input_file=\$INPUT

# Create a temporary output directory specific to this job
temp_output="$temp_path/\${sample}_haplotype_\${SLURM_JOB_ID}"
mkdir -p "\$temp_output"

# Define the path to the reference genome
GENOME="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/unzipped/GRCh38.primary.genome.fa"

# Define the path for the output GVCF file
OUTPUT="\$temp_output/\$(basename "\$input_file" _sorted.dup.bqsr.bam).g.vcf.gz"

# Specify Docker container images for GATK and Samtools
GATK_IMAGE="/mnt/beegfs/apptainer/images/gatk_latest.sif"
GATK_IMAGE="/mnt/beegfs/apptainer/images/gatk4.sif"
SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"

# Print important variables for debugging
echo "Input file: \$input_file"
echo "Sample: \$sample"
echo "Output: \$OUTPUT"

# Check the integrity of the input BAM file
srun apptainer exec "\$SAMTOOLS_IMAGE" samtools quickcheck "\$input_file"

if [ \$? -eq 0 ]; then
    echo "The BAM file \$input_file is OK."
else
    echo "The BAM file \$input_file failed the quick check and may be corrupted."
fi

#check if @RG SM is correct now since i needed to correct before in ReagGroup_change >_<
srun apptainer exec "\$SAMTOOLS_IMAGE" samtools view -H \$input_file | grep '@RG'

# check if index of bam exists 
if [ ! -f "\$input_file.bai" ]; then
    srun apptainer exec "\$SAMTOOLS_IMAGE" samtools index -@ 40 "\$input_file"
else
    echo "\$input_file.bai already exists."
fi

# Run HaplotypeCaller within the GATK Docker container
srun apptainer exec "\$GATK_IMAGE" gatk --java-options "-Xmx120G -XX:ParallelGCThreads=40 -XX:ConcGCThreads=40" HaplotypeCaller \
    -I "\$input_file" \
    -R "\$GENOME" \
    -ERC GVCF \
    -O "\$OUTPUT" \
    --tmp-dir "\$temp_output"

# Create the final output directory and move the results
if [ \$? -eq 0 ]; then
    echo "HaplotypeCaller completed successfully."

    # Create the final output directory and move the results
    output_dest_DIR=\$(dirname "\$INPUT" | sed 's/bam\$/VarCal/')
    mkdir -p "\$output_dest_DIR"
    mv "\$OUTPUT" "\$output_dest_DIR/"
    mv "\${OUTPUT}.tbi" "\$output_dest_DIR/"
else
    echo "Error: HaplotypeCaller failed."
fi






EOL