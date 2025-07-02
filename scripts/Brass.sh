normal_bam=(
    /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/align_newref/C2040neg/bam/C2040neg_sorted.dup.bqsr.bam
            /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/align_newref/C2041neg/bam/C2041neg_sorted.dup.bqsr.bam
            )

tumor_bam=(
    /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/align_newref/C2040ttd/bam/C2040ttd_sorted.dup.bqsr.bam
            /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/align_newref/C2041ttd/bam/C2041ttd_sorted.dup.bqsr.bam
            )

all_bam=(
    /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/align_newref/C2040neg/bam/C2040neg_sorted.dup.bqsr.bam
    /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/align_newref/C2041neg/bam/C2041neg_sorted.dup.bqsr.bam
    /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/align_newref/C2040ttd/bam/C2040ttd_sorted.dup.bqsr.bam
    /mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/align_newref/C2041ttd/bam/C2041ttd_sorted.dup.bqsr.bam
    )

sbatch_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/scripts/brass
mkdir -p $sbatch_path

reference_fa="/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/reference/genome/unzipped/GRCh38.primary.genome.fa"


temp_path="/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/temp/brass"
log_path=/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/output/logs/brass
mkdir -p $log_path


mkdir -p "$temp_path"

brass_image=/mnt/beegfs/apptainer/images/dockstore-cgpwgs_latest.sif

echo "Creating scripts..."
setup_sbatch="$sbatch_path/01_bam_stats.sbatch"

cat >"$setup_sbatch" <<EOL
#!/bin/bash
#SBATCH --job-name=bam_stats
#SBATCH --output=$log_path/01_brass_%A_%a.out
#SBATCH --time=72:00:00
#SBATCH --array=0-3%4
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=40


brass=$brass_image

all_bam=(${all_bam[@]})


input="\${all_bam[\$SLURM_ARRAY_TASK_ID]}"

output_path=\${input}.bas

echo "creating bam stats for \$input"
apptainer exec \$brass bam_stats -i \$input -o \$output_path -@ 36


EOL


setup_sbatch="$sbatch_path/02_bamtoBigwig.sbatch"

cat >"$setup_sbatch" <<EOL
#!/bin/bash
#SBATCH --job-name=bamtoBw
#SBATCH --output=$log_path/02_bamtoBw_%A_%a.out
#SBATCH --time=72:00:00
#SBATCH --array=0-3%4
#SBATCH --nodes=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40


if [[ \$input_normal == *"C2040"* ]]; then
    sample="C2040"
else
    sample="C2041"
fi

reference_genome=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/reference/\${sample}neg_reference.fa


brass="$brass_image"

all_bam=(${all_bam[@]})

input="\${all_bam[\$SLURM_ARRAY_TASK_ID]}"

output_path="\${input}.bw"

sample_name=\$(basename "\$input" .bam).bw

temp_path="/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/temp/brass/bigwig"

mkdir -p "\$temp_path"

temp_output="\$temp_path/\$sample_name"

echo "Creating bam stats for \$input"
echo "temp output \$temp_output"
##apptainer exec "\$brass" bamToBw.pl -t 36 -b "\$input" -o "\$temp_output" -r \$reference_genome
#apptainer exec "\$brass" bam2bw -i "\$input" -o "\$temp_output" -r \$reference_genome

final_dir=\$(dirname \$(dirname "\$input"))/bigWig

mkdir -p "\$final_dir"

#mv "\$temp_output" "\$final_dir"

bigwig=\$final_dir/\$sample_name
output_dir=\$final_dir/depth
mkdir -p \$output_dir

echo "Finding high coverage positions" 
apptainer exec "\$brass" detectExtremeDepth -b \$bigwig -o \$output_dir


# Directory previous to input
#final_dir=\$(dirname \$(dirname "\$input"))/bigWig

#mkdir -p "\$final_dir"

#mv "\$temp_output" "\$final_dir"

EOL

echo "Script created at: $setup_sbatch"




###########################################

setup_sbatch="$sbatch_path/03_brass.sbatch"

cat >"$setup_sbatch" <<EOL
#!/bin/bash
#SBATCH --job-name=brass
#SBATCH --output=$log_path/03_brass_%A_%a.out
#SBATCH --time=72:00:00
#SBATCH --array=0-1%4
#SBATCH --nodes=1
#SBATCH --mem=240G
#SBATCH --cpus-per-task=40



brass=$brass_image

normal=(${normal_bam[@]})
tumor=(${tumor_bam[@]})


input_normal="\${normal[\$SLURM_ARRAY_TASK_ID]}"
input_tumor="\${tumor[\$SLURM_ARRAY_TASK_ID]}"

if [[ \$input_normal == *"C2040"* ]]; then
    sample="C2040"
else
    sample="C2041"
fi

temp_path=${temp_path}/brass_\${sample}
reference_genome=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/reference/\${sample}neg_reference.fa
depth=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/align_newref/\${sample}neg/bigWig/depth/\${sample}neg_sorted.dup.bqsr.bed

temp_path="/mnt/beegfs/scratch/MCFONSECA/pedro.silvestre/brca_project/temp/brass"

output_path=\$(dirname \$input_normal)


echo "creating bam stats for \$input"

apptainer exec \$brass brass.pl \\
    -tumour \$input_tumor \\
    -normal \$input_normal \\
    -depth \$depth \\
    -genome \$reference_genome \\
    -species human \\
    -protocol WGS \\
    -outdir \$temp_path


#    brass.pl [options]

#      Required parameters:
#        -outdir    -o   Folder to output result to.
#        -tumour    -t   Tumour BAM file
#        -normal    -n   Normal BAM file
#        -depth     -d   Regions of excessive sequencing depth to be ignored
#       -genome    -g   Genome fasta file
#        -species   -s   Species name
#        -assembly  -as  Assembly name
#        -protocol  -pr  Sequencing protocol (WGS|WXS|RNA)
#        -g_cache   -gc  Genome cache file.
#        -viral     -vi  Virus sequences from NCBI
#        -microbe   -mi  Microbe sequences file prefix from NCBI, please exclude '.N.fa.2bit'
#        -gcbins    -b   5 column bed coord file, col 4 number of non-N bases, col 5 GC fraction.
#        -cytoband    -cb  Cytoband file for a species build (can be obtained from UCSC)
#        -centtel   -ct  TSV file of usable regions of the chromosome
#                          Example in perl/share/Rdefault/

#      Optional
#        -mingroup  -j   Minimum reads to call group [2].
#        -mincn       -cn   Minimum CN change for copynumber_flag [0.3].
#        -repeats   -r   Repeat file, see 'make-repeat-file' (legacy)
#        -sampstat  -ss  ASCAT sample statistics file or file containing
#                          rho 0.XXXXX [0.75] (~ 1-normalContamination)
#                          Ploidy X.XXX [2.0]
#                          GenderChr Y [Y]
#                          GenderChrFound Y/N [Y]
#        -platform    -pl  Sequencing platform (when not defined in BAM header)

#        -tum_name    -tn  Tumour sample name (when not defined in BAM header)
#        -norm_name   -nn  Normal sample name (when not defined in BAM header)
#        -filter      -f   bgzip tabix-ed normal panel groups file
#        -noclean     -x   Don't tidyup the processing areas.
#        -cpus        -c   Number of cores to use. [1]
#                         - recommend max 2 during 'input' process.

#      Targeted processing (further detail under PROCESS OPTIONS):
#        -process   -p   Only process this step then exit, optionally set -index
#        -index     -i   Optionally restrict '-p' to single job
#        -limit     -l   Define with -p and -i, see below


final_dir=/mnt/nfs/lobo/MCFONSECA-NFS/pedro.silvestre/BRCA_project/align_newref/StructVariants_brass
mkdir -p \$final_dir

#mv \$temp_path \$final_dir

EOL
