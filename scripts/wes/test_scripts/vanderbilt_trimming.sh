#!/bin/bash
#SBATCH --job-name=fastp_trim_vanderbilt
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/fastp_%j.out
#SBATCH --error=logs/fastp_%j.err

module purge
module load bear-apps/2023a
module load fastp/0.24.0-GCC-12.3.0

RAW_DIR="/rds/projects/e/elhamsak-group5/main_project/vanderbilt_samples"
OUT_DIR="/rds/projects/e/elhamsak-group5/main_project/trimmed_files/vanderbilt_samples"

for SAMPLE in $(cat /rds/projects/e/elhamsak-group5/main_project/vanderbilt_samples/paper_samples.txt); do
  R1="${RAW_DIR}/${SAMPLE}_1.fastq.gz"
  R2="${RAW_DIR}/${SAMPLE}_2.fastq.gz"

  echo "▶ Processing $SAMPLE"

  fastp \
    -i "$R1" \
    -I "$R2" \
    -o "${OUT_DIR}/${SAMPLE}_1.trimmed.fastq.gz" \
    -O "${OUT_DIR}/${SAMPLE}_2.trimmed.fastq.gz" \
    -h "${OUT_DIR}/${SAMPLE}_fastp.html" \
    -j "${OUT_DIR}/${SAMPLE}_fastp.json" \
    -w 8
done

