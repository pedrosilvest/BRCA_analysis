#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples
script_name=FilterVCF_control
script_number=01
mem=240G
cpus_per_task=40


# list of realpaths for files

example=(
   /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/vcf/by_sample/C2041neg.vcf.gz
   /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/vcf/by_sample/C2040neg.vcf.gz
)

#change this
files_path=(
   /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples/vcf/2.by_sample/CBE7negB.vcf.gz
   )

total_arrays=$(( ${#files_path[@]} - 1 ))

# FILTER CONDITIONS: 

# Define common criteria for most variants

common_criteria="(FMT/DP >= 30) & (FMT/GQ > 30)" 

# Criteria for heterozygous variants with AD[1] > 10 -> if previous failed, then check if alt allele AD is higher than 13 FMT/AD[sample:allele] check alt allele 1 for the first sample 0

heterozygous_criteria="((FMT/GT = \"1/0\" | FMT/GT = \"0/1\" | FMT/GT = \"0|1\" | FMT/GT= \"1|0\") & MIN(FMT/AD) >= 20)" # baixar  mexer neste

# Criteria for homozygous alternate variants -> if it is hom ref and one of the alt AD is lower then 7.
hom_alt_criteria="((GT == \"1/1\" | GT == \"1|1\") & FMT/AD[0:1] >= 30)"


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

sample=\$(basename "\$INPUT" .vcf.gz)

temp_output="$temp_path/\${sample}_filterVCF_\${SLURM_JOB_ID}"
mkdir -p "\$temp_output"

tmp_file="\$temp_output/\$(basename "\$INPUT" .vcf.gz)_tmp.vcf.gz"

final_file="\$temp_output/\$(basename "\$INPUT" .vcf.gz)_for_fasta.vcf.gz"

######### IMAGE #########
BCFTOOLS=$BCFTOOLS
########

GENOME=$GENOME

echo "\$INPUT"
echo "\$tmp_file"
echo "\$final_file"

# FILTER CONDITIONS: 

# Define common criteria for most variants

common_criteria="(FMT/DP >= 30) & (FMT/GQ > 30)" 

# Criteria for heterozygous variants with AD[1] > 10 -> if previous failed, then check if alt allele AD is higher than 13 FMT/AD[sample:allele] check alt allele 1 for the first sample 0

heterozygous_criteria="((FMT/GT = \"1/0\" | FMT/GT = \"0/1\" | FMT/GT = \"0|1\" | FMT/GT= \"1|0\") & MIN(FMT/AD) >= 20)" # baixar  mexer neste

# Criteria for homozygous alternate variants -> if it is hom ref and one of the alt AD is lower then 7.
#old
hom_alt_criteria="((GT == \"1/1\" | GT == \"1|1\") & FMT/AD[0:1] >= 30)"


srun apptainer exec \$BCFTOOLS bcftools norm -f \$GENOME -m-any \$INPUT --threads 30 -Oz -o \$tmp_file

srun apptainer exec \$BCFTOOLS bcftools filter -i "\$final_criteria" --threads 30 \$tmp_file -Oz -o \$final_file


### TESTING RATIOS

final_path=$output_path/vcf/3.for_fasta
mkdir -p \$final_path

TABIX=$TABIX

srun apptainer exec "\$TABIX" tabix -p vcf "\$final_file"

source /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/miniconda3/etc/profile.d/conda.sh

plot_ratios=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/python_script/obtain_ratio_variant/plot_ratios.py

apptainer exec "\$BCFTOOLS" bcftools query -l "\$final_file" | while read -r sample; do
    if [[ "\$final_file" == *no_hSOR* ]]; then
        prefix="no_hSOR_"
    else
        prefix=""
    fi

    ratio_file="\$final_path/\${prefix}\${sample}_ratio_all_variants.txt"
    png_file="\$final_path/\${prefix}\${sample}_ratio_variants.png"

    apptainer exec "\$BCFTOOLS" bcftools query -i 'GT!="0/0" && GT!="0|0" && GT!=".|." && GT!="./."' -s "\$sample" -f '[%SAMPLE\t%GT\t%AD\t%DP\n]' "\$final_file" | awk -F'\t' '{
        split(\$3, AD, ",");
        ratio = AD[2] / \$4;
        print \$0 "\t" ratio;
    }' > "\$ratio_file"

    python "\$plot_ratios" --input "\$ratio_file" --output "\$png_file"

done

#INDEX VCF


srun apptainer exec "\$TABIX" tabix -p vcf "\$final_file"


### MOVING
mv \$final_file* \$final_path

rm -rf \$temp_output


EOL