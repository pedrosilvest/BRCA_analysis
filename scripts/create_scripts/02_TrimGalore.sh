#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/CBE7
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project
script_name=TrimGalore
script_number=02
mem=240G
cpus_per_task=40


# Should be a list of paths!!
fasta_paths=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/file_path_CBE7.txt


#create folders
logs_path=$scratch_path/output/logs
scripts_path=$scratch_path/scripts
temp_path=$scratch_path/temp

#mkdir -p $output_path
mkdir -p $logs_path
mkdir -p $scripts_path
mkdir -p $temp_path


# arrays setup
total_arrays=$(( $(cat "$fasta_paths" | wc -l) - 1 ))



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

file_path=$fasta_paths

IFS=$'\r\n' GLOBIGNORE='*' command eval 'INPUTS=(\$(cat "\$file_path"))'

INPUT="\${INPUTS[\$SLURM_ARRAY_TASK_ID]}"

# Extract sample name from the file path
sample=\$(basename "\$INPUT" | cut -d "_" -f 1)

echo \$sample


# Create a directory for the sample in a temporary folder with job ID
temp_output="$temp_path/\${sample}_TrimGalore_\${SLURM_JOB_ID}"
mkdir -p "\$temp_output"

#R1 and R2 files


# image
FASTQC="/mnt/beegfs/apptainer/images/mcfonsecalab_fastqc_preprocessing_latest.sif"
over_represented_sequence=GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
#over_represented_sequence=CGGTGGTCGCCGTATCATTAGATCGGAAGAGCACACGTCTGAACTCCAGT

# Replace R1 with R2 in the file name
R1=\$INPUT
R2=\$(echo "\$INPUT" | sed 's/1.fq/2.fq/')
echo "\$R1"
# Run TrimGalore on both R1 and R2 in the temporary folder 
#srun apptainer exec "\$FASTQC" fastqc -t 20 -o "\$temp_output" "\$R1" "\$R2"
srun apptainer exec "\$FASTQC" trim_galore --cores 8 --paired --adapter \$over_represented_sequence --length 70 --clip_R1 1 --clip_R2 1 --quality 30 --fastqc -o \$temp_output \$R1 \$R2
#removes adapter
# Move the contents of the temporary folder to the final destination

sample_dest_folder="$output_path/\$sample/trim"

mkdir -p \$sample_dest_folder

mv "\$temp_output"/* \$sample_dest_folder

echo "final dest: \$sample_dest_folder"


# Remove the temporary folder
rmdir "\$temp_output"


EOL