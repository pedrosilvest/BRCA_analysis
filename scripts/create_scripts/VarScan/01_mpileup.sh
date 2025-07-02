#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/VarScan
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples
script_name=MpileUp
script_number=01
mem=240G
cpus_per_task=40

# IMAGES
BWA="/mnt/beegfs/apptainer/images/bwa_latest.sif"
SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"

GENOME="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/unzipped/GRCh38.primary.genome.fa"


#create folders
logs_path=$scratch_path/output/logs
scripts_path=$scratch_path/scripts
temp_path=$scratch_path/temp

#mkdir -p $output_path
mkdir -p $logs_path
mkdir -p $scripts_path
mkdir -p $temp_path


total_arrays=$(find "$output_path"/*/bam -type f -name "*dup.bqsr.bam" ! -path "*CBE7*" | wc -l)
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

GENOME=$GENOME

file_path=\$(find "$output_path"/*/bam -type f -name "*dup.bqsr.bam" ! -path "*CBE7*")


IFS=$'\r\n' GLOBIGNORE='*' command eval 'INPUTS=(\$(echo "\$file_path"))'

INPUT="\${INPUTS[\$SLURM_ARRAY_TASK_ID]}"

# bam file
bam_file="\$INPUT"



sample=\$(basename "\$INPUT" | cut -d "_" -f 1)


# Create a directory for the sample in a temporary folder with job ID
temp_output="$temp_path/\${sample}_mpileup_\${SLURM_JOB_ID}"
mkdir -p "\$temp_output"



# IMAGES
SAMTOOLS_IMAGE=$SAMTOOLS_IMAGE

pileup_file=\$temp_output/\${sample}.pileup

srun apptainer exec \$SAMTOOLS_IMAGE samtools mpileup \\
    -f $GENOME \\
    -o \$pileup_file \\
    \$bam_file


mpile_dir=$output_path/\$sample/mpileup

mkdir -p \$mpile_dir




mv \$pileup_file \$mpile_dir

rm -rf "\$temp_output"



EOL