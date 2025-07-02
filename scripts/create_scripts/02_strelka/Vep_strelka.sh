#!/bin/bash
#SBATCH --job-name=variants_gene_count
#SBATCH --output=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/out_%A_%a.out
#SBATCH --array=0
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40

# Define Variables
VEP_IMAGE="/mnt/beegfs/apptainer/images/ensembl-vep_release_108.2.sif"
CACHE_DIR="/mnt/nfs/lobo/IMM-NFS/ensembl_vep/cache/"
INPUT_DIR="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/strelka"
OUTPUT_DIR="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/strelka"
SAMPLES=("C2040" "C2041")

# Function to process each sample
run_vep() {
    local sample=$1
    local input_file="${INPUT_DIR}/${sample}/StrelkaSomaticWorkflow/results/variants/somatic.snvs.vcf.gz"
    local output_file="${OUTPUT_DIR}/${sample}/StrelkaSomaticWorkflow/results/variants/vep_somatic.snvs.vcf.gz"

    if [[ ! -f "$input_file" ]]; then
        echo "Error: Input file $input_file not found!" >&2
        return 1
    fi

    echo "Processing sample: $sample"
    echo "Input: $input_file"
    echo "Output: $output_file"

    start_time=$(date +%s)
    
    apptainer exec "$VEP_IMAGE" vep \
        --cache --dir "$CACHE_DIR" \
        --canonical --vcf --force_overwrite \
        --hgvs --af_gnomadg -i "$input_file" -o "$output_file"

    end_time=$(date +%s)
    duration=$((end_time - start_time))
    echo "VEP processing for $sample completed in $duration seconds."
}

# Process each sample
for sample in "${SAMPLES[@]}"; do
    run_vep "$sample"
done

echo "All samples processed successfully."
