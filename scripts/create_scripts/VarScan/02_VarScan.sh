#Setup paths
scratch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/VarScan
output_path=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/samples
script_name=VarScan
script_number=02
mem=240G
cpus_per_task=40

# IMAGES
varscan_image=/mnt/beegfs/apptainer/images/varscan.sif

GENOME="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/unzipped/GRCh38.primary.genome.fa"

samples=("C2040" "C2041") # "CBE7")

#create folders
logs_path=$scratch_path/output/logs
scripts_path=$scratch_path/scripts
temp_path=$scratch_path/temp

#mkdir -p $output_path
mkdir -p $logs_path
mkdir -p $scripts_path
mkdir -p $temp_path


total_arrays=${#samples[@]} 
total_arrays=$((total_arrays - 1))


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

GENOME=$GENOME


INPUTS=(${samples[@]})
sample="\${INPUTS[\$SLURM_ARRAY_TASK_ID]}"

file_path=\$(find "$output_path"/*/mpileup -type f -name "\${sample}*pileup")

normal_pileup=""
tumor_pileup=""

# Assign files based on their names
for file in \$file_path; do
    if [[ "\$file" == *"neg"* ]]; then
        normal_pileup="\$file"
    elif [[ "\$file" == *"ttd"* ]]; then
        tumor_pileup="\$file"
    fi
done



# Create a directory for the sample in a temporary folder with job ID
temp_output="$temp_path/\${sample}_mpileup_\${SLURM_JOB_ID}"
mkdir -p "\$temp_output"



srun apptainer exec $varscan_image java "-Xmx120G" "-XX:ParallelGCThreads=36" "-XX:ConcGCThreads=36" -jar /opt/varscan/VarScan.jar somatic \\
    \$normal_pileup \\
    \$tumor_pileup \\
    \$temp_output/\$sample

#VarScan.jar somatic normal.pileup tumor.pileup output.basename


varscan_dir=$output_path/VarScan

mkdir -p \$varscan_dir




#mv \$temp_output* \$varscan_dir

#rm -rf "\$temp_output"



EOL





#USAGE: VarScan somatic [normal_pileup] [tumor_pileup] [Opt: output] OPTIONS
#	normal_pileup - The SAMtools pileup file for Normal
#	tumor_pileup - The SAMtools pileup file for Tumor
#	output - Output base name for SNP and indel output

#OPTIONS:
#	--output-snp - Output file for SNP calls [output.snp]
#	--output-indel - Output file for indel calls [output.indel]
#	--min-coverage - Minimum coverage in normal and tumor to call variant [8]
#	--min-coverage-normal - Minimum coverage in normal to call somatic [8]
#	--min-coverage-tumor - Minimum coverage in tumor to call somatic [6]
#	--min-var-freq - Minimum variant frequency to call a heterozygote [0.10]
#	--min-freq-for-hom	Minimum frequency to call homozygote [0.75]
#	--normal-purity - Estimated purity (non-tumor content) of normal sample [1.00]
#	--tumor-purity - Estimated purity (tumor content) of tumor sample [1.00]
#	--p-value - P-value threshold to call a heterozygote [0.99]
#	--somatic-p-value - P-value threshold to call a somatic site [0.05]
#	--strand-filter - If set to 1, removes variants with >90% strand bias [0]
#	--validation - If set to 1, outputs all compared positions even if non-variant
#	--output-vcf - If set to 1, output VCF instead of VarScan native format


# Step 3: Post-process SNPs and Indels
#java -jar VarScan.jar processSomatic output.snp
#java -jar VarScan.jar processSomatic output.indel