#!/bin/bash

# Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples
script_name=Manta
script_number=00
mem=240G
cpus_per_task=40

normal_bam="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples/CBE7negB/bam/CBE7negB_sorted.dup.bqsr.bam"
tumoral_bam="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples/CBE7ttdB/bam/CBE7ttdB_sorted.dup.bqsr.bam"
genome_ref=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/unzipped/GRCh38.primary.genome.fa

# Create folders
logs_path=$scratch_path/output/logs/03_Manta
scripts_path=$scratch_path/scripts/03_Manta
temp_path=$scratch_path/temp

mkdir -p $logs_path
mkdir -p $scripts_path
mkdir -p $temp_path

setup_sbatch="$scripts_path/${script_number}_${script_name}.sbatch"
echo "Creating scripts... $setup_sbatch"

manta_image=/mnt/beegfs/apptainer/images/manta_v1.6.0.sif

cat > "$setup_sbatch" <<EOL
#!/bin/bash
#SBATCH --job-name=manta
#SBATCH --output=$logs_path/01_manta_%A_%a.out
#SBATCH --time=72:00:00
#SBATCH --array=0
#SBATCH --nodes=1
#SBATCH --mem=$mem
#SBATCH --cpus-per-task=$cpus_per_task

manta_image=$manta_image
normal_bam="$normal_bam"
tumoral_bam="$tumoral_bam"

timestamp() {
    date +"%Y-%m-%d %H:%M:%S"
}

export -f timestamp

temp_output=$temp_path/CBE7

echo "\$(timestamp) -> Job started! Configuring a new workflow running script.."
#srun apptainer exec \$manta_image configManta.py \\
#    --normalBam \$normal_bam \\
#    --tumorBam \$tumoral_bam \\
#    --referenceFasta $genome_ref \\
#    --runDir \$temp_output

srun apptainer exec \$manta_image configManta.py \\
    --normalBam \$normal_bam \\
    --tumorBam \$tumoral_bam \\
    --referenceFasta $genome_ref \\
    --runDir \$temp_output \\
    --minEdgeObservations 2 \\
    --minCandidateVariantScore 5 \\
    --minCandidateSpanningCount 2 \\
    --minCandidateClusterSize 2 \\
    --enable-duplicate-edges




echo "\$(timestamp) -> Configuration completed. Will run now workflow."
srun apptainer exec \$manta_image \$temp_output/runWorkflow.py -j 30

echo -e "\$(timestamp) -> Finished job."
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
echo -e "\$(timestamp) -> All done!"

# Fixed incorrect variable reference for moving results
final_path=$output_path/manta_SV_lower
mkdir -p \$final_path

mv \$temp_output/* \$final_path/
EOL
