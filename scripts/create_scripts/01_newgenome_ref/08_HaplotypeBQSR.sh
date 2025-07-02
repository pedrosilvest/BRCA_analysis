#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples/newref
script_name=HaplotypeBQSR
script_number=08
mem=240G
cpus_per_task=40


# IMAGE
GATK_IMAGE="/mnt/beegfs/apptainer/images/gatk_latest.sif"
SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"
GENOME=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples/reference/CBE7negB_reference.fa


#create folders
logs_path=$scratch_path/output/logs/01_newgenome_ref
scripts_path=$scratch_path/scripts/01_newgenome_ref
temp_path=$scratch_path/temp

#mkdir -p $output_path
mkdir -p $logs_path
mkdir -p $scripts_path
mkdir -p $temp_path


total_arrays=$(find $output_path/*/bam -type f -name "*sorted.bam" | wc -l)
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


# this will run the haplotype caller to produce a vcf instead of gvcf that will later be used to produce a vcf of known variants

file_path=\$(find $output_path/*/bam -type f -name "*_sorted.dup.bam")

echo "all file paths: \$file_path"

IFS=$'\r\n' GLOBIGNORE='*' command eval 'INPUTS=(\$(echo "\$file_path"))'


INPUT="\${INPUTS[\$SLURM_ARRAY_TASK_ID]}"


sample=\$(basename "\$INPUT" _sorted.dup.bam)


input_file=\$INPUT
echo "Input file: \$input_file"
echo "Sample: \$sample"

echo "\$input_file"

GATK_IMAGE=$GATK_IMAGE
SAMTOOLS_IMAGE=$SAMTOOLS_IMAGE

GENOME=$GENOME


# Check if index of BAM exists, if not, create it
if [ ! -f "\$input_file.bai" ]; then
    srun apptainer exec "\$SAMTOOLS_IMAGE" samtools index -@ 40 "\$input_file"
else
    echo "\$input_file.bai already exists."
fi



# Create temporary output directory
temp_output="$temp_path/\${sample}_haplotypeBQSR_\${SLURM_JOB_ID}"
mkdir -p "\$temp_output"

# Define output paths

OUTPUT="\$temp_output/\$(basename "\$input_file" _sorted.dup.bam)_raw.vcf.gz"

# Run HaplotypeCaller within GATK Docker container
srun apptainer exec "\$GATK_IMAGE" gatk --java-options "-Xmx200G -XX:ParallelGCThreads=35 -XX:ConcGCThreads=35" HaplotypeCaller \
    -I "\$input_file" \
    -R "\$GENOME" \
    -O "\$OUTPUT" \
    --tmp-dir "\$temp_output"
    #-ERC GVCF

# Move output files to appropriate directories


final_dir="$output_path/KnowVar"

echo "moving to \$final_dir"

mkdir -p "\$final_dir"

mv \$OUTPUT* "\$final_dir/"

rm -rf \$temp_output

EOL
