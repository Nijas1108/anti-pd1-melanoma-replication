#!/bin/bash
#SBATCH --job-name=star_align
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=logs/star_align_%j.out
#SBATCH --error=logs/star_align_%j.err

set -euo pipefail

mkdir -p logs

cd /rds/projects/e/elhamsak-group5/main_project/

module purge
module load bear-apps/2024a
module load STAR/2.7.11b-GCC-13.3.0

# Paths
GENOME_DIR="/rds/projects/e/elhamsak-group5/main_project/reference_genome/STAR_index"
INPUT_DIR="/rds/projects/e/elhamsak-group5/main_project/processed_data/trimmed_files/rnaseq"
OUT_DIR="/rds/projects/e/elhamsak-group5/main_project/processed_data/aligned_files/rnaseq"


# Your paired-end SRR list
SRR_LIST=(
"SRR3184279"
"SRR3184280"
"SRR3184281"
"SRR3184282"
"SRR3184283"
"SRR3184284"
"SRR3184285"
"SRR3184286"
"SRR3184287"
"SRR3184288"
"SRR3184289"
"SRR3184290"
"SRR3184291"
"SRR3184292"
"SRR3184293"
"SRR3184294"
"SRR3184295"
"SRR3184296"
"SRR3184297"
"SRR3184298"
"SRR3184299"
"SRR3184300"
"SRR3184301"
"SRR3184302"
"SRR3184303"
"SRR3184304"
"SRR3184305"
"SRR3184306"
)

for SRR in "${SRR_LIST[@]}"; do
  echo "========================================"
  echo "Aligning $SRR with STAR"
  echo "========================================"

  R1="${INPUT_DIR}/${SRR}_1.trimmed.fastq.gz"
  R2="${INPUT_DIR}/${SRR}_2.trimmed.fastq.gz"
  PREFIX="${OUT_DIR}/${SRR}_"

  if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "ERROR: Missing trimmed FASTQs for $SRR"
    echo "  $R1"
    echo "  $R2"
    exit 1
  fi

  STAR \
    --runThreadN "${SLURM_CPUS_PER_TASK}" \
    --genomeDir "$GENOME_DIR" \
    --readFilesIn "$R1" "$R2" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$PREFIX" \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --twopassMode Basic \
    --outSAMattrRGline ID:"$SRR" SM:"$SRR" PL:ILLUMINA \
    --outSAMstrandField intronMotif

  echo "Done: $SRR"
done

echo "All STAR alignments completed."