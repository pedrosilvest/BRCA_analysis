#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples
script_name=ChangeGT
script_number=02
mem=240G
cpus_per_task=40


# list of realpaths for files

example=(
   /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/vcf/by_sample/C2041neg.vcf.gz
   /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/vcf/by_sample/C2040neg.vcf.gz
)

#change this
files_path=(
   /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples/vcf/3.for_fasta/CBE7negB_for_fasta.vcf.gz
   )

total_arrays=$(( ${#files_path[@]} - 1 ))

# IMAGE
GATK_IMAGE="/mnt/beegfs/apptainer/images/gatk4.sif"
SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"
TABIX="/mnt/beegfs/apptainer/images/mcfonsecalab_variantutils_latest.sif"
BCFTOOLS=/mnt/beegfs/apptainer/images/mcfonsecalab_variantutils_latest.sif

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

sample=\$(basename "\$INPUT" _for_fasta.vcf.gz)

temp_output="$output_path/\${sample}_newvcf_\${SLURM_JOB_ID}"
mkdir -p "\$temp_output"


output_file="\$temp_output/\${sample}_noHet_for_fasta.vcf.gz"
final_file="\$temp_output/\${sample}_noHetnoDel_for_fasta.vcf.gz"

source /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/miniconda3/etc/profile.d/conda.sh

conda activate

script_python=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/python_script/processVCF/VCFchangeGT.py

python \$script_python -i \$INPUT -o \$output_file


######### IMAGE #########
BCFTOOLS=$BCFTOOLS
########

apptainer exec \$BCFTOOLS bcftools view -i 'ALT!="*"' \$output_file -o \$final_file


TABIX=$TABIX

srun apptainer exec "\$TABIX" tabix -p vcf "\$output_file"
srun apptainer exec "\$TABIX" tabix -p vcf "\$final_file"


final_path=$output_path/vcf/3.for_fasta
mkdir -p \$final_path

### MOVING
mv \$final_file* \$final_path
mv \$output_file* \$final_path

rm -rf \$temp_path


EOL

