#!/bin/bash
#SBATCH --job-name=mutect2_ucla
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=1-1
#SBATCH --output=logs/mutect2_%A_%a.out
#SBATCH --error=logs/mutect2_%A_%a.err

set -euo pipefail

# Load modules
module purge
module load bluebear
module load bear-apps/2022b
module load GATK/4.4.0.0-GCCcore-12.2.0-Java-17

#-------------------------------
# Define your samples as a Bash array
#-------------------------------
SAMPLES=(
SRR4289727
)

#-------------------------------
# Get sample for this SLURM array task
#-------------------------------
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}
BAM=${SAMPLE}.sorted.bam

#-------------------------------
# Paths
#-------------------------------
REF=/rds/projects/e/elhamsak-group5/main_project/reference_genome/hg38.fa
TUMOR_BAM=/rds/projects/e/elhamsak-group5/main_project/aligned_files/UCLA_samples/${SAMPLE}.sorted.bam

OUT_DIR=/rds/projects/e/elhamsak-group5/main_project/variant_calling/UCLA_samples/${SAMPLE}
mkdir -p "$OUT_DIR" logs

REF=/rds/projects/e/elhamsak-group5/main_project/reference_genome/hg38.fa
GNOMAD=/path/to/af-only-gnomad.vcf.gz
COMMON_SNP=/path/to/common_biallelic_snps.vcf.gz
# Optional for WES (recommended). If you don’t have it, set TARGETS="" and remove -L lines below.
TARGETS=/path/to/targets.interval_list

# Basic input checks (fail fast, helpful on clusters)
[[ -f "${BAM}" ]] || { echo "Missing BAM: ${BAM}"; exit 1; }
[[ -f "${BAM}.bai" || -f "${BAM}.sorted.bam.bai" || -f "${SAMPLE}.sorted.bam.bai" ]] || true

# ---- Mutect2 (tumor-only) ----
gatk Mutect2 \
  -R "${REF}" \
  -I "${BAM}" \
  -tumor "${SAMPLE}" \
  -L "${TARGETS}" \
  --germline-resource "${GNOMAD}" \
  -O "vcf/${SAMPLE}.unfiltered.vcf.gz"

# ---- Contamination estimation ----
gatk GetPileupSummaries \
  -I "${BAM}" \
  -V "${COMMON_SNP}" \
  -L "${TARGETS}" \
  -O "tables/${SAMPLE}.pileups.table"

gatk CalculateContamination \
  -I "tables/${SAMPLE}.pileups.table" \
  -O "tables/${SAMPLE}.contamination.table" \
  --tumor-segmentation "tables/${SAMPLE}.segments.table"

# ---- Filter calls ----
gatk FilterMutectCalls \
  -R "${REF}" \
  -V "vcf/${SAMPLE}.unfiltered.vcf.gz" \
  --contamination-table "tables/${SAMPLE}.contamination.table" \
  --tumor-segmentation "tables/${SAMPLE}.segments.table" \
  -O "vcf/${SAMPLE}.filtered.vcf.gz"

echo "[$(date)] Done: ${SAMPLE}"
