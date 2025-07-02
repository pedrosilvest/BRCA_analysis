#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples
script_name=GenotypeGVCFs
script_number=11
mem=240G
cpus_per_task=40

######### IMAGE #########
GATK_IMAGE="/mnt/beegfs/apptainer/images/gatk_latest.sif"
SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"
########

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
#SBATCH --array=0-24%4
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --mem=$mem
#SBATCH --cpus-per-task=$cpus_per_task

######### IMAGE #########
GATK_IMAGE=$GATK_IMAGE
SAMTOOLS_IMAGE=$SAMTOOLS_IMAGE
########


# Function to determine the chromosome name based on task ID
get_chromosome_name() {
    case "\$1" in
        23) echo "chrX" ;;
        24) echo "chrY" ;;
        25) echo "chrM" ;;
        *) echo "chr\$((\$1))" ;;
    esac
}

# Determine the chromosome name
CHROMOSOME=\$(get_chromosome_name \$((SLURM_ARRAY_TASK_ID + 1)))

# Define your paths



WORKSPACE_PATH="gendb://$output_path/genomicsdb/genomicsdb_\$CHROMOSOME"
GENOME_REFERENCE=$GENOME


temp_output="$temp_path/genomicsdb_\${CHROMOSOME}"

# Create the output directory if it doesn't exist
mkdir -p "\$temp_output"


echo \$WORKSPACE_PATH

# Generate a per-chromosome output VCF file
OUTPUT_Chromosome="\$temp_output/\${CHROMOSOME}.vcf"

# Run GenotypeGVCFs for the specified chromosome
srun apptainer exec "\$GATK_IMAGE" gatk --java-options "-Xmx200G -XX:ParallelGCThreads=40 -XX:ConcGCThreads=40" GenotypeGVCFs \
    -R \$GENOME_REFERENCE \
    --tmp-dir \$temp_output \
    -V \$WORKSPACE_PATH \
    -O \$OUTPUT_Chromosome


final_dir=$output_path/vcf

mkdir -p \$final_dir
mkdir -p \$final_dir/0.by_chr

mv \$OUTPUT_Chromosome \$final_dir/0.by_chr

rm -rf \$temp_output

EOL