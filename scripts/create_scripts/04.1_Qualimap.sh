#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples
script_name=qualimap
script_number=04.1
mem=240G
cpus_per_task=40

# IMAGES
samtools_image="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"


#create folders
logs_path=$scratch_path/output/logs
scripts_path=$scratch_path/scripts
temp_path=$scratch_path/temp

#mkdir -p $output_path
mkdir -p $logs_path
mkdir -p $scripts_path
mkdir -p $temp_path


total_arrays=$(find $output_path/*/bam -type f -name "*.bam" ! -name "*.dup.bqsr.bam" | wc -l)
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


file_path=\$(find $output_path/*/bam -type f -name "*.bam" ! -name "*.dup.bqsr.bam")

samtools_image=$samtools_image
IFS=$'\r\n' GLOBIGNORE='*' command eval 'INPUTS=(\$(echo "\$file_path"))'

INPUT="\${INPUTS[\$SLURM_ARRAY_TASK_ID]}"

sample=\$(basename "\$INPUT" .bam)


# Create a directory for the sample in a temporary folder with job ID
temp_output="$temp_path/\${sample}_qualimap_\${SLURM_JOB_ID}"
mkdir -p "\$temp_output"


echo "sample name: \$sample"

image_path=/mnt/beegfs/apptainer/images/qualimap_2.2.1.sif


srun apptainer exec \$image_path qualimap bamqc \\
        -bam \$INPUT \\
        -outdir \$temp_output \\
        -gd HUMAN \\
        -nt \$SLURM_CPUS_PER_TASK \\
        --java-mem-size=24G


mkdir -p $output_path/\${sample}/qualimap


mv "\$temp_output"/* $output_path/\${sample}/qualimap

rm -rf "\$temp_output"

EOL


