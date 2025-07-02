#!/bin/bash

# Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples
script_name=Gridss
script_number=00
mem=240G
cpus_per_task=40

normal_bam="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples/CBE7negB/bam/CBE7negB_sorted.dup.bqsr.bam"
tumoral_bam="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples/CBE7ttdB/bam/CBE7ttdB_sorted.dup.bqsr.bam"
genome_ref=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/unzipped/GRCh38.primary.genome.fa

# Create folders
logs_path=$scratch_path/output/logs/04_Gridss
scripts_path=$scratch_path/scripts/04_Gridss
temp_path=$scratch_path/temp

mkdir -p $logs_path
mkdir -p $scripts_path
mkdir -p $temp_path

setup_sbatch="$scripts_path/${script_number}_${script_name}.sbatch"
echo "Creating scripts... $setup_sbatch"


grids_image=/mnt/beegfs/apptainer/images/gridss_latest.sif

cat > "$setup_sbatch" <<EOL
#!/bin/bash
#SBATCH --job-name=$script_name
#SBATCH --output=$logs_path/01_gridss_%A_%a.out
#SBATCH --time=72:00:00
#SBATCH --array=0
#SBATCH --nodes=1
#SBATCH --mem=$mem
#SBATCH --cpus-per-task=$cpus_per_task

grids_image=$grids_image

temp_output=$temp_path/CBE7

srun apptainer exec $grids_image gridss \\
  -r $genome_ref \\
  -j /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/scripts/create_scripts/04_gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar \\
  -o \$temp_output/all_calls.vcf \\
  -t 30 \\
  -w \$temp_output \\
  $normal_bam \\
  $tumoral_bam

srun apptainer exec $grids_image gridss_somatic_filter \
  --pondir refdata/hg19/dbs/gridss/pon3792v1/ \
  --input \$temp_output/all_calls.vcf \
  --output high_confidence_somatic.vcf.gz \
  --fulloutput high_and_low_confidence_somatic.vcf.gz \
  --scriptdir $(dirname $(which gridss_somatic_filter)) \
  -n 1 \
  -t 2


EOL
