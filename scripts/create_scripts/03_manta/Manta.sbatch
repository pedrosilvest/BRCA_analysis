#!/bin/bash
#SBATCH --job-name=manta
#SBATCH --output=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7/temp/manta_%A_%a.out
#SBATCH --time=72:00:00
#SBATCH --array=0
#SBATCH --nodes=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40


normal_bam="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples/CBE7negB/bam/CBE7negB_sorted.dup.bqsr.bam"
tumoral_bam="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples/CBE7ttdB/bam/CBE7ttdB_sorted.dup.bqsr.bam"
genome_ref=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/unzipped/GRCh38.primary.genome.fa

manta_image=/mnt/beegfs/apptainer/images/manta_v1.6.0.sif
export -f timestamp

echo "$(timestamp) -> Job started! Configuring a new workflow running script.."
srun apptainer exec /mnt/beegfs/apptainer/images/manta_v1.6.0.sif configManta.py --normalBam $normal_bam --tumorBam $tumoral_bam --referenceFasta $genome_ref --runDir /mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/temp/mantaC2040

echo "$(timestamp) -> Configuration completed. Will run now workflow."
# -j is cpus
srun apptainer exec $manta_image /mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7/temp/mantaCBE7/runWorkflow.py -j 30

echo -e "$(timestamp) -> Finished job."
echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID
echo -e "$(timestamp) -> All done!"



mkdir -p /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples/manta
mv /mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/temp/mantaC2040 /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples/manta


