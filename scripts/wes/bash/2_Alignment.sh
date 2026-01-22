#!/bin/bash
#SBATCH --job-name=alignment
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/bwa_%A_%a.out
#SBATCH --error=logs/bwa_%A_%a.err
#SBATCH --array=0-15 # adjust array to sample number - 1

set -euo pipefail

## decide what samples you want to run this for vanderbilt or ucla
#default vanderbilt

# Load modules
module purge
module load bluebear
module load bear-apps/2023a
module load BWA/0.7.17-GCCcore-12.3.0
module load SAMtools/1.21-GCC-12.3.0

#General paths
REF_FASTA="/rds/projects/e/elhamsak-group5/main_project/reference_genome/hg38/hg38.fa"
SCRATCH_DIR=/scratch/$USER

# Paths UCLA
#FASTQ_DIR="/rds/projects/e/elhamsak-group5/main_project/processed_data/trimmed_files/UCLA_samples"
#OUT_DIR="/rds/projects/e/elhamsak-group5/main_project/processed_data/aligned_files/UCLA_samples"


# Paths Vanderbilt
FASTQ_DIR="/rds/projects/e/elhamsak-group5/main_project/processed_data/trimmed_files/vanderbilt_samples"
OUT_DIR="/rds/projects/e/elhamsak-group5/main_project/processed_data/aligned_files/vanderbilt_samples"


mkdir -p logs "$OUT_DIR" "$SCRATCH_DIR"

# List of UCLA sample names
SRR_LIST=(
#SRR3083837 SRR3083838 SRR3083839 SRR3083840 SRR3083841 SRR3083842
#SRR3083843 SRR3083844 SRR3083845 SRR3083846 SRR3083847 SRR3083848
#SRR3083849 SRR3083850 SRR3083851 SRR3083852 SRR3083853 SRR3083854
#SRR3083855 SRR3083856 SRR3083857 SRR3083858 SRR3083859 SRR3083860
#SRR3083861 SRR3083862 SRR3083863 SRR3083864 SRR3083866 SRR3083867
#SRR3083868 SRR3083869 SRR3083870 SRR3083871 SRR3083872 SRR3083873
#SRR3083874 SRR3083875 SRR3083876 SRR3083877 SRR3083878 SRR3083879
#SRR3083880 SRR3083881 SRR3083882 SRR3083883 SRR4289714 SRR4289715
#SRR4289716 SRR4289717 SRR4289718 SRR4289719 SRR4289720 SRR4289721
#SRR4289722 SRR4289723 SRR4289724 SRR4289725 SRR4289726 
)

#List of Vanderbilt sample names
SRR_LIST=(
SRR4289728 SRR4289730 SRR4289732 SRR4289734 SRR4289736 SRR4289738
SRR4289740 SRR4289744 SRR4289727 SRR4289729 SRR4289731 SRR4289733
SRR4289735 SRR4289737 SRR4289739 SRR4289743
)

# Get the sample for this job array index
SAMPLE=${SRR_LIST[$SLURM_ARRAY_TASK_ID]}
R1="${FASTQ_DIR}/${SAMPLE}_1.trimmed.fastq.gz"
R2="${FASTQ_DIR}/${SAMPLE}_2.trimmed.fastq.gz"
SORTED_BAM="${OUT_DIR}/${SAMPLE}.sorted.bam"
TMP_PREFIX="${SCRATCH_DIR}/${SAMPLE}.sort.tmp"

# Skip if already done
if [[ -f "$SORTED_BAM" ]]; then
    echo "Sorted BAM already exists for $SAMPLE. Skipping."
    exit 0
fi

# Input checks
if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "ERROR: Missing FASTQ files for $SAMPLE"
    exit 1
fi

echo "Aligning sample: $SAMPLE"
echo "Using $SLURM_CPUS_PER_TASK CPUs"
echo "Temporary files in: $TMP_PREFIX"

# Alignment + sorting (optimised)
# -R adds header 
bwa mem -t $SLURM_CPUS_PER_TASK \
  -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:lib1\tPL:ILLUMINA\tPU:${SAMPLE}.lane1" \
  "$REF_FASTA" "$R1" "$R2" | \
samtools view -b - | \
# -l - compression level (up to 9)
samtools sort \
  -@ $SLURM_CPUS_PER_TASK \
  -m 3G \
  -l 1 \
  -T "$TMP_PREFIX" \
  -o "$SORTED_BAM"

# Index BAM
samtools index "$SORTED_BAM"

echo "Alignment and indexing completed for $SAMPLE"
