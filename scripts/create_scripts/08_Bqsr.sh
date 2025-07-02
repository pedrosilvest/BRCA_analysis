#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project
script_name=Bqsr
script_number=08
mem=240G
cpus_per_task=40



GENOME="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/GRCh38.primary.genome.fa.gz"

#create folders
logs_path=$scratch_path/output/logs
scripts_path=$scratch_path/scripts
temp_path=$scratch_path/temp

#mkdir -p $output_path
mkdir -p $logs_path
mkdir -p $scripts_path
mkdir -p $temp_path


file_paths=(
"/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/CBE7negB"
"/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/CBE7ttdB"
)

# Get the total number of files (array length)
total_files=${#file_paths[@]}

# Subtract 1 from the total number of files
total_arrays=$((total_files - 1))


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


file_paths=($(printf '"%s"\n' "${file_paths[@]}"))

INPUT="\${file_paths[\$SLURM_ARRAY_TASK_ID]}"
echo "\$INPUT"

sample=\$(basename "\$INPUT")

INPUT=(\$INPUT/bam/\${sample}_sorted.dup.bam)

input_file=\$INPUT
echo "Input file: \$input_file"
echo "Sample: \$sample"


temp_output="$temp_path/\${sample}_\${SLURM_JOB_ID}_bqsr"
mkdir -p "\$temp_output"


#GENOME AND KNOWN LOCATION
GENOME="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/unzipped/GRCh38.primary.genome.fa"
KNOWN="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/known-sites/Homo_sapiens_assembly38.dbsnp138.vcf"

#BASE RECAL OUTPUT TABLE
TABLE=\$temp_output/\$(basename "\$input_file" .bam).recal_data.table

#TEMPORARY FILES
TMP=\$temp_output

#BQSR OUTPUT
OUTPUT=\$temp_output/\$(basename "\$input_file" .bam).bqsr.bam


echo "input: \$input_file"
echo "output table: \$TABLE"
echo "output: \$OUTPUT"



######### IMAGE #########
GATK_IMAGE="/mnt/beegfs/apptainer/images/gatk_latest.sif"
SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"
########


#CHECK INPUT TO COPY
srun apptainer exec "\$SAMTOOLS_IMAGE" samtools quickcheck "\$input_file"

if [ \$? -eq 0 ]; then
    echo "The BAM file \$input_file is ok"
else
    echo "The BAM file \$input_file failed the quick check and may be corrupted."
fi


START_TIME=\$(date +%s)
echo "BASE RECAL: \$input_file"

srun apptainer exec "\$GATK_IMAGE" gatk --java-options "-Xmx120G -XX:ParallelGCThreads=40 -XX:ConcGCThreads=40" BaseRecalibrator \
    -I "\$input_file" \
    -R "\$GENOME" \
    --known-sites "\$KNOWN" \
    -O "\$TABLE" \
    --tmp-dir "\$TMP"

echo "BQSR: \$input_file"
echo "bqsr table: \$TABLE"
srun apptainer exec "\$GATK_IMAGE" gatk --java-options "-Xmx120G -XX:ParallelGCThreads=40 -XX:ConcGCThreads=40" ApplyBQSR \
    -I "\$input_file" \
    -R "\$GENOME" \
    --bqsr-recal-file "\$TABLE" \
    -O "\$OUTPUT" \
    --tmp-dir "\$TMP"


srun apptainer exec "\$SAMTOOLS_IMAGE" samtools quickcheck "\$OUTPUT"

srun apptainer exec "\$SAMTOOLS_IMAGE" samtools index -@ 40 "\$OUTPUT"


output_bam_DIR=\$(dirname "\$INPUT")

bqsr_table_DIR=\$(dirname "\$INPUT" | sed 's/bam\$/bqsr_table/')
mkdir -p \$bqsr_table_DIR


echo "final paths: bam - \$output_bam_DIR table - \$bqsr_table_DIR"

if [ \$? -eq 0 ]; then
    echo "EVERYTHING IS FINE"
    mv \$TABLE \$bqsr_table_DIR
    mv \$temp_output/\$sample* \$output_bam_DIR
else
    echo "The BAM file \$OUTPUT failed the quick check and may be corrupted."
fi





EOL


