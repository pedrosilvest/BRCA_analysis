#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples
script_name=GenomicsDB
script_number=10.1
mem=240G
cpus_per_task=40


# number of arrays = 24 since it is the number of chr including Y

g_vcf_files=$(find "$output_path"/*/VarCal -type f -name "*.g.vcf.gz" | sed 's/^/\t-V /' | sed 's/$/ \\/')




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
#SBATCH --array=0-24%4
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --mem=$mem
#SBATCH --cpus-per-task=$cpus_per_task


######### IMAGE #########
GATK_IMAGE="/mnt/beegfs/apptainer/images/gatk_latest.sif"
SAMTOOLS_IMAGE="/mnt/beegfs/apptainer/images/mcfonsecalab_htstools_plus_latest.sif"
########

GENOME="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/unzipped/GRCh38.primary.genome.fa"

# Check if the chromosome is 23 (for human sex chromosomes X)
get_chromosome_name() {
    case "\$1" in
        23) echo "chrX" ;;
        24) echo "chrY" ;;
        25) echo "chrM" ;;
        *) echo "chr\$((\$1))" ;;
    esac
}

CHROMOSOME=\$(get_chromosome_name \$((SLURM_ARRAY_TASK_ID + 1)))

mkdir -p $output_path/genomicsdb/

WORKSPACE_PATH="$output_path/genomicsdb/genomicsdb_\$CHROMOSOME/"

temp_output="/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/temp/genomicsdb_\${CHROMOSOME}"

mkdir -p "\$temp_output"

echo \$CHROMOSOME
echo \$WORKSPACE_PATH

srun apptainer exec "\$GATK_IMAGE" gatk --java-options "-Xmx120G -XX:ParallelGCThreads=40 -XX:ConcGCThreads=40" GenomicsDBImport \\
  -R /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/unzipped/GRCh38.primary.genome.fa \\
  --genomicsdb-workspace-path \$WORKSPACE_PATH \\
  --tmp-dir \$temp_output \\
  $g_vcf_files
  -L \$CHROMOSOME




EOL

