#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/featureCounts_%j.out
#SBATCH --error=logs/featureCounts_%j.err

set -euo pipefail

mkdir -p logs

module purge
module load bear-apps/2022a
module load Subread/2.0.4-GCC-11.3.0

GTF="/rds/projects/e/elhamsak-group5/main_project/reference_genome/STAR_index/gencode.v49.annotation.nochr.fixed.gtf"
BAM_DIR="/rds/projects/e/elhamsak-group5/main_project/processed_data/aligned_files/rnaseq"
OUT_DIR="/rds/projects/e/elhamsak-group5/main_project/processed_data/rna_seq/rnaseq_counts"
OUT="${OUT_DIR}/gene_counts.txt"

featureCounts \
  -T "${SLURM_CPUS_PER_TASK}" \
  -s 0 \
  -p --countReadPairs \
  -B -C \
  -t exon \
  -g gene_id \
  -a "$GTF" \
  -o "$OUT" \
  "$BAM_DIR"/*_Aligned.sortedByCoord.out.bam