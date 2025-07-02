normal_bam=(/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/C2040neg/bam/C2040neg_sorted.dup.bqsr.bam
            /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/C2041neg/bam/C2041neg_sorted.dup.bqsr.bam)

tumor_bam=(/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/C2040ttd/bam/C2040ttd_sorted.dup.bqsr.bam
            /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/C2041ttd/bam/C2041ttd_sorted.dup.bqsr.bam)

sbatch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/scripts/manta
mkdir -p $sbatch_path

final_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/Manta
mkdir -p $final_path

for i in "${!normal_bam[@]}"; do


    reference_genome="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/unzipped/GRCh38.primary.genome.fa"

    temp_path="/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/temp/mantaC204${i}"
    log_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/output/logs/manta/C204${i}
    mkdir -p $log_path
    mkdir -p "$temp_path"

    manta_image=/mnt/beegfs/apptainer/images/manta_v1.6.0.sif

    echo "Creating scripts..."
    setup_sbatch="$sbatch_path/01_manta_C204${i}_setup.sbatch"

    cat >"$setup_sbatch" <<EOL
#!/bin/bash
#SBATCH --job-name=manta
#SBATCH --output=$log_path/01_manta_%A_%a.out
#SBATCH --time=72:00:00
#SBATCH --array=0
#SBATCH --nodes=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40


normal_bam="${normal_bam[i]}"
tumoral_bam="${tumor_bam[i]}"

manta_image=$manta_image
export -f timestamp

echo "\$(timestamp) -> Job started! Configuring a new workflow running script.."
srun apptainer exec $manta_image configManta.py \
--normalBam \$normal_bam \
--tumorBam \$tumoral_bam \
--referenceFasta $reference_genome \
--runDir $temp_path

echo "\$(timestamp) -> Configuration completed. Will run now workflow."
# -j is cpus
srun apptainer exec \$manta_image $temp_path/runWorkflow.py -j 30

echo -e "\$(timestamp) -> Finished job."
echo "Statistics for job \$SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j \$SLURM_JOB_ID
echo -e "\$(timestamp) -> All done!"




mv $temp_path $final_path


EOL

sbatch $setup_sbatch
done