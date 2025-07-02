#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples
script_name=Align
script_number=04
mem=240G
cpus_per_task=40

# IMAGES
BWA="/mnt/beegfs/apptainer/images/bwa_latest.sif"
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


total_arrays=$(find $output_path/*/trim -type f -name "*_1_*fq.gz" | wc -l)
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

file_path=\$(find $output_path/*/trim -type f -name "*1.fq.gz")


IFS=$'\r\n' GLOBIGNORE='*' command eval 'INPUTS=(\$(echo "\$file_path"))'

INPUT="\${INPUTS[\$SLURM_ARRAY_TASK_ID]}"

# R1 and R2 files
R1="\$INPUT"
R2="\${R1//_1_/_2_}"
R2="\${R2//val_1/val_2}"

echo "R1: \$R1"
echo "R2: \$R2"

sample=\$(basename "\$INPUT" | cut -d "_" -f 1)
full_name=\$(basename "\$R1" | awk -F '_1_' '{print \$1}')


# Create a directory for the sample in a temporary folder with job ID
temp_output="$temp_path/\${sample}_align_\${SLURM_JOB_ID}"
mkdir -p "\$temp_output"



BAM_OUTPUT="\$temp_output/\$full_name.bam"  # Updated this line
SAVE="\$temp_output/\$full_name""_sorted.bam"
STATS_OUTPUT="\$temp_output/\$full_name""_stats.txt"

echo "full_name: \$full_name"
echo "folder name: \$sample"
echo "sorted_output: \$SAVE"
echo "bam output: \$BAM_OUTPUT"
echo "stats output: \$STATS_OUTPUT"

SM="\$sample"
ID="\$full_name""_id"
RG="@RG\tID:\$ID\tSM:\$SM\tPL:ILLUMINA"    # AQUI DEVIA TER TIDO CUIDADO POR TER LANES. O SM é a SAMPLE C2040neg e o ID é o unique id.

echo "\$RG"

# IMAGES
BWA=$BWA
SAMTOOLS_IMAGE=$SAMTOOLS_IMAGE

# BAM MEM ALIGNMENT

apptainer exec "\$BWA" bwa mem -t 30 -R "\$RG" "\$GENOME" "\$R1" "\$R2" | \
  apptainer exec "\$SAMTOOLS_IMAGE" samtools view -b -o "\$BAM_OUTPUT" -

apptainer exec "\$SAMTOOLS_IMAGE" samtools sort -@ 30 \$BAM_OUTPUT -o "\$SAVE"

# STATISTICS
apptainer exec "\$SAMTOOLS_IMAGE" samtools flagstat -@ 30 "\$SAVE" > "\$STATS_OUTPUT"

apptainer exec "\$SAMTOOLS_IMAGE" samtools view -H "\$SAVE" | grep "^@RG"

mkdir -p $output_path/newref
mkdir -p $output_path/newref/\$sample

bam_dir="$output_path/newref/\$sample/bam/"
stats_dir="$output_path/newref/\$sample/stats/"
mkdir -p \$bam_dir
mkdir -p \$stats_dir

mv "\$temp_output"/*sorted.bam "\$bam_dir/"
mv "\$temp_output"/*stats.txt "\$stats_dir"

rm -rf "\$temp_output"



EOL