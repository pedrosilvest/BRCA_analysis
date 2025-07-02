#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples/newref
script_name=Bqsr
script_number=09
mem=240G
cpus_per_task=40

# arrays setup
total_arrays=$(find $output_path/*/bam -type f -name "*sorted.dup.bam" | wc -l)
total_arrays=$((total_arrays - 1))

GENOME=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples/reference/CBE7negB_reference.fa
KNOWN="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples/vcf/3.for_fasta/CBE7negB_for_fasta.vcf.gz"



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

file_path=\$(find $output_path/*/bam -type f -name "*_sorted.dup.bam")

echo "all file paths: \$file_path"

IFS=$'\r\n' GLOBIGNORE='*' command eval 'INPUTS=(\$(echo "\$file_path"))'


INPUT="\${INPUTS[\$SLURM_ARRAY_TASK_ID]}"

sample=\$(basename "\$INPUT" _sorted.dup.bam)


input_file=\$INPUT

echo "Input file: \$input_file"
echo "Sample: \$sample"



temp_output="$temp_path/\${sample}_\${SLURM_JOB_ID}_bqsr"
mkdir -p "\$temp_output"


#GENOME AND KNOWN LOCATION
GENOME=$GENOME
KNOWN=$KNOWN

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


