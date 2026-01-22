#!/bin/bash
#SBATCH --job-name=fastp_trim
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/fastp_%j.out
#SBATCH --error=logs/fastp_%j.err

module purge
module load bear-apps/2023a
module load fastp/0.24.0-GCC-12.3.0

INPUT_DIR="/rds/projects/e/elhamsak-group5/main_project/rnaseq"
OUT_DIR="/rds/projects/e/elhamsak-group5/main_project/processed_data/trimmed_files/rnaseq"


SRR_LIST=(
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
  echo "=============================="
  echo "Running fastp for $SRR"
  echo "=============================="

  R1="${INPUT_DIR}/${SRR}_1.fastq.gz"
  R2="${INPUT_DIR}/${SRR}_2.fastq.gz"

  OUT_R1="${OUT_DIR}/${SRR}_1.trimmed.fastq.gz"
  OUT_R2="${OUT_DIR}/${SRR}_2.trimmed.fastq.gz"

  HTML="${OUT_DIR}/${SRR}.fastp.html"
  JSON="${OUT_DIR}/${SRR}.fastp.json"

  fastp \
    -i "$R1" -I "$R2" \
    -o "$OUT_R1" -O "$OUT_R2" \
    -h "$HTML" -j "$JSON" \
    -w 8

  echo "Done: $SRR"
done

echo "All fastp jobs completed."