#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples
script_name=Merge_chrVCF
script_number=12
mem=240G
cpus_per_task=40

######### IMAGE #########
GATK_IMAGE="/mnt/beegfs/apptainer/images/gatk4.sif"
SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"
########

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
#SBATCH --array=0
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --mem=$mem
#SBATCH --cpus-per-task=$cpus_per_task


input_directory="$output_path/vcf/0.by_chr"

final_dir="$output_path/vcf/0.by_chr"

temp_output="$temp_path/mergeVCF_chr_\${SLURM_JOB_ID}"

# Create the output directory if it doesn't exist
mkdir -p "\$temp_output"

find \$input_directory -name "chr*.vcf" > \$temp_output/input_vcf_files.list


##############################

GATK_IMAGE=$GATK_IMAGE
SAMTOOLS_IMAGE=$SAMTOOLS_IMAGE

##############################



srun apptainer exec "\$GATK_IMAGE" gatk MergeVcfs \
    -I \$temp_output/input_vcf_files.list \
    -O \$temp_output/all_chr.vcf

mv \$temp_output/all_chr.vcf \$final_dir

rm -rf \$temp_output


EOL
