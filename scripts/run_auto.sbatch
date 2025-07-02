#!/bin/bash
#SBATCH --job-name=runAuto
#SBATCH --output=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/scripts/postjob_%j.out
#SBATCH --time=01:00:00


# Your actual job commands go here
echo "Running post-processing job after Job IDs"

chmod +x /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/scripts/create_scripts/01_newgenome_ref/06_FilterBam.sh

./mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/scripts/create_scripts/01_newgenome_ref/06_FilterBam.sh

sbatch /mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7/scripts/01_newgenome_ref/06_FilterBam.sh