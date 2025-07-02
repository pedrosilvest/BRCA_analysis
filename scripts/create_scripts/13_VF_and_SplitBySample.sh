#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples
script_name=VF_SplitBySample
script_number=13
mem=240G
cpus_per_task=40

# IMAGE
GATK_IMAGE="/mnt/beegfs/apptainer/images/gatk4.sif"
SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"
TABIX="/mnt/beegfs/apptainer/images/mcfonsecalab_variantutils_latest.sif"
BCFTOOLS=/mnt/beegfs/apptainer/images/mcfonsecalab_variantutils_latest.sif

GENOME="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/unzipped/GRCh38.primary.genome.fa"

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



# ----------------- VARIANT FILTRATION PIPELINE ------------------------------ #

echo "Variant Filtration Pipeline starting..."
##################
GATK_IMAGE=$GATK_IMAGE
SAMTOOLS_IMAGE=$SAMTOOLS_IMAGE
TABIX=$TABIX
BCFTOOLS=$BCFTOOLS
##################


### TABIX FILE

GENOME=$GENOME

INPUT=$output_path/vcf/0.by_chr/all_chr.vcf

VF_output_path=$output_path/vcf/1.VariantFiltration


srun apptainer exec "\$TABIX" tabix -p vcf "\$INPUT"

temp_output=$temp_path/VariantFiltration_\${SLURM_JOB_ID}

# Create the output directory if it doesn't exist
mkdir -p \$temp_output

name=\$(basename "\$INPUT" .vcf.gz)

OUTFILE_FILTERED="\$temp_output/\${name}_vf_filtered.vcf"

OUTFILE_SELECTED="\$temp_output/\${name}_sv_filtered.vcf"

srun apptainer exec "\$GATK_IMAGE" gatk --java-options "-Xmx120G -XX:ParallelGCThreads=40 -XX:ConcGCThreads=40" VariantFiltration \
   -R "\$GENOME" \
   -V "\$INPUT" \
   --filter-name LowQD --filter-expression "QD < 2.0" \
   --filter-name LowMQ --filter-expression "MQ < 40.0" \
   --filter-name HighFS --filter-expression "FS > 60.0" \
   --filter-name LowMQRankSum --filter-expression "MQRankSum < -12.5" \
   --filter-name HighSOR --filter-expression "SOR > 3.0" \
   --filter-name LowReadPosRankSum --filter-expression "ReadPosRankSum < -8.0" \
   -O "\$OUTFILE_FILTERED"

srun apptainer exec "\$GATK_IMAGE" gatk --java-options "-Xmx120G -XX:ParallelGCThreads=40 -XX:ConcGCThreads=40" SelectVariants \
   -R "\$GENOME" \
   -V "\$OUTFILE_FILTERED" \
   --exclude-non-variants \
   --exclude-filtered \
   -O "\$OUTFILE_SELECTED"


mkdir -p \$VF_output_path

mv \$OUTFILE_FILTERED \$VF_output_path
mv \$OUTFILE_SELECTED \$VF_output_path

rm -rf \$temp_output



# ------------------ Split By Sample Pipeline ------------------------ #

Split_output_path=$output_path/vcf/2.by_sample

echo "SplitBySample Pipeline starting..."


INPUT="\$VF_output_path/\${name}_sv_filtered.vcf"

# Capture the sample names in an array
SAMPLES=(\$(srun apptainer exec \$BCFTOOLS bcftools query -l "\$INPUT"))


temp_output=$temp_path/SplitVCF_sample_\${SLURM_JOB_ID}
# Create the output directory if it doesn't exist
mkdir -p \$temp_output

mkdir -p "\$Split_output_path"

# split the VCF file
for SAMPLE_NAME in "\${SAMPLES[@]}"; do
  echo "Processing sample: \$SAMPLE_NAME"
  OUTPUT="\$temp_output/\$SAMPLE_NAME.vcf.gz"
  srun apptainer exec \$BCFTOOLS bcftools view -O z -s "\$SAMPLE_NAME" "\$INPUT" -o "\$OUTPUT"
  mv \$OUTPUT \$Split_output_path
  echo "Finished processing sample: \$SAMPLE_NAME"
done

EOL