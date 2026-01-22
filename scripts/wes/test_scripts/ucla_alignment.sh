#!/bin/bash
#SBATCH --job-name=alignment_ucla
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/bwa_%A_%a.out
#SBATCH --error=logs/bwa_%A_%a.err
#SBATCH --array=0-0

set -euo pipefail

# Load modules
module purge
module load bluebear
module load bear-apps/2023a
module load BWA/0.7.17-GCCcore-12.3.0
module load SAMtools/1.21-GCC-12.3.0

# Paths
FASTQ_DIR=/rds/projects/e/elhamsak-group5/main_project/trimmed_files/UCLA_samples
REF_FASTA=/rds/projects/e/elhamsak-group5/main_project/reference_genome/hg38.fa
OUT_DIR=/rds/projects/e/elhamsak-group5/main_project/aligned_files/UCLA_samples
SCRATCH_DIR=/scratch/$USER

mkdir -p logs "$OUT_DIR" "$SCRATCH_DIR"

# List of sample names
SRR_LIST=(
SRR4289724
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
bwa mem -t $SLURM_CPUS_PER_TASK \
  -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:lib1\tPL:ILLUMINA\tPU:${SAMPLE}.lane1" \
  "$REF_FASTA" "$R1" "$R2" | \
samtools view -b - | \
samtools sort \
  -@ $SLURM_CPUS_PER_TASK \
  -m 3G \
  -l 1 \
  -T "$TMP_PREFIX" \
  -o "$SORTED_BAM"

# Index BAM
samtools index "$SORTED_BAM"

echo "Alignment and indexing completed for $SAMPLE"