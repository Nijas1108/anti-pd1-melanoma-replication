#!/bin/bash
#SBATCH --job-name=alignment_vanderbilt
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/bwa_%j.out
#SBATCH --error=logs/bwa_%j.err

# Load modules
module purge
module load bluebear
module load bear-apps/2023a
module load BWA/0.7.17-GCCcore-12.3.0
module load SAMtools/1.21-GCC-12.3.0

# Paths
FASTQ_DIR=/rds/projects/e/elhamsak-group5/main_project/trimmed_files/vanderbilt_samples
REF=/rds/projects/e/elhamsak-group5/main_project/reference_genome/hg38.fa
OUT_DIR=/rds/projects/e/elhamsak-group5/main_project/aligned_files/vanderbilt_samples

mkdir -p logs
mkdir -p "$OUT_DIR"

THREADS="${SLURM_CPUS_PER_TASK}"
SORT_THREADS=2

for R1 in "$FASTQ_DIR"/*_1.trimmed.fastq.gz; do
  SAMPLE=$(basename "$R1" _1.trimmed.fastq.gz)
  R2="$FASTQ_DIR/${SAMPLE}_2.trimmed.fastq.gz"
  BAM="$OUT_DIR/${SAMPLE}.sorted.bam"

  if [[ ! -f "$R2" ]]; then
    echo "Skipping $SAMPLE missing R2"
    continue
  fi

  echo "Aligning sample: $SAMPLE"
#bwa + defining read groups
  bwa mem -t "$THREADS" \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:lib1\tPL:ILLUMINA\tPU:${SAMPLE}.lane1" \
    "$REF" "$R1" "$R2" \
    | samtools sort -@ "$SORT_THREADS" -o "$BAM"

  samtools index "$BAM"
done