#!/bin/bash
grids_image=/mnt/beegfs/apptainer/images/gridss_latest.sif



normal_bam=(/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/C2040neg/bam/C2040neg_sorted.dup.bqsr.bam
            /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/C2041neg/bam/C2041neg_sorted.dup.bqsr.bam)

tumor_bam=(/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/C2040ttd/bam/C2040ttd_sorted.dup.bqsr.bam
            /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/C2041ttd/bam/C2041ttd_sorted.dup.bqsr.bam)

sbatch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/scripts/grids
mkdir -p $sbatch_path


for i in "${!normal_bam[@]}"; do


    reference_genome="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/unzipped/GRCh38.primary.genome.fa"

    temp_path="/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/temp/gridssC204${i}"
    log_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/output/logs/gridss/C204${i}
    mkdir -p $log_path
    mkdir -p "$temp_path"

    manta_image=/mnt/beegfs/apptainer/images/manta_v1.6.0.sif

    grids_image=/mnt/beegfs/apptainer/images/gridss_latest.sif

    echo "Creating scripts..."
    setup_sbatch="$sbatch_path/01_gridss_C204${i}_setup.sbatch"

    cat >"$setup_sbatch" <<EOL
#!/bin/bash
#SBATCH --job-name=grids
#SBATCH --output=$log_path/01_manta_%A_%a.out
#SBATCH --time=72:00:00
#SBATCH --array=0
#SBATCH --nodes=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40


normal_bam="${normal_bam[i]}"
tumoral_bam="${tumor_bam[i]}"



srun apptainer exec $grids_image gridss \\
  -r $reference_genome \\
  -j /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/scripts/create_scripts/04_gridss/gridss-2.13.2-gridss-jar-with-dependencies.jar \\
  -o \$temp_output/all_calls.vcf \\
  -t 30 \\
  -w $temp_path \\
  \$normal_bam \\
  \$tumoral_bam



EOL

sbatch $setup_sbatch
done