#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples
script_name=NewFasta
script_number=03
mem=240G
cpus_per_task=40


# list of realpaths for files

example=(
   /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/vcf/by_sample/C2041neg.vcf.gz
   /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/vcf/by_sample/C2040neg.vcf.gz
)

#change this
files_path=(
   /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples/vcf/3.for_fasta/CBE7negB_noHetnoDel_for_fasta.vcf.gz
   )

total_arrays=$(( ${#files_path[@]} - 1 ))

# IMAGE
GATK_IMAGE="/mnt/beegfs/apptainer/images/gatk4.sif"
SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"
TABIX="/mnt/beegfs/apptainer/images/mcfonsecalab_variantutils_latest.sif"
BCFTOOLS=/mnt/beegfs/apptainer/images/mcfonsecalab_variantutils_latest.sif
BWA="/mnt/beegfs/apptainer/images/bwa_latest.sif"

GENOME="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/unzipped/GRCh38.primary.genome.fa"

#create folders
logs_path=$scratch_path/output/logs/01_newgenome_ref
scripts_path=$scratch_path/scripts/01_newgenome_ref
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

# Store the file paths in an array

files_path="${files_path[@]}"

IFS=$' ' GLOBIGNORE='*' command eval 'INPUTS=(\$files_path)'

# output will be SAMPLE_NAME_filtered.vcf.gz

# Get the specific input file based on the SLURM_ARRAY_TASK_ID
INPUT="\${INPUTS[\$SLURM_ARRAY_TASK_ID]}"

sample=\$(basename "\$INPUT" _noHetnoDel_for_fasta.vcf.gz)

temp_output="$temp_path/\${sample}_newfasta"
mkdir -p "\$temp_output"

chain_file="\$temp_output/\${sample}_reference.chain"

final_file="\$temp_output/\${sample}_reference.fa"

######### IMAGE #########
BCFTOOLS=$BCFTOOLS
########


GENOME=$GENOME

#echo "creating new fasta:  \$pre_final"

apptainer exec "$BCFTOOLS" bcftools consensus -f "\$GENOME" \$INPUT --chain \$chain_file -o \$final_file




## --------------------------------------------

### CREATE INDEX AND DICTIONARY FOR NEW FASTA

SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"

GATK_IMAGE=$GATK_IMAGE

BWA=$BWA

# Set the output directory
output_dir="$output_path/reference"

mkdir -p \$output_dir
mkdir -p \$output_dir/liftover

# Index the genome file with samtools
srun apptainer exec \$SAMTOOLS_IMAGE samtools faidx \$final_file

# Create a sequence dictionary

srun apptainer exec "\$GATK_IMAGE" gatk CreateSequenceDictionary R="\$final_file" O="\${final_file%.fa}.dict"

srun apptainer exec \$BWA bwa index "\$final_file"

echo "Genome processing completed successfully for \$final_file."


## ------- Move 
mv \$temp_output/*.chain \$output_dir/liftover
mv \$temp_output/* \$output_dir





EOL
