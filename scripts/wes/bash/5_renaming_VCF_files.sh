#!/bin/bash
#SBATCH --job-name=rename_chromosomes
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/rename_chromosomes_%j.out
#SBATCH --error=logs/rename_chromosomes_%j.err

set -euo pipefail

module purge
module load bluebear
module load bear-apps/2024a
module load BCFtools/1.21-GCC-13.3.0

# -------------------------------
# Inputs PATHS
# -------------------------------
#General
chr_map_file="/rds/projects/e/elhamsak-group5/main_project/reference_genome/nochr_to_chr.txt"

#Vanderbilt
BASE_DIR="/rds/projects/e/elhamsak-group5/main_project/processed_data/wes/variant_calling/vanderbilt_processed_bed"

#UCLA
#BASE_DIR="/rds/projects/e/elhamsak-group5/main_project/processed_data/wes/variant_calling/ucla_processed_bed"


# Check if the chromosome mapping file exists
if [[ ! -f "$chr_map_file" ]]; then
  echo "Chromosome mapping file not found: $chr_map_file"
  exit 1
fi

echo "Renaming chromosomes in VCF files under directory: $BASE_DIR"
echo "Using chromosome mapping file: $chr_map_file"

# Loop through all subdirectories under BASE_DIR that contain filtered.PASS.vcf.gz files
for VCF_DIR in "$BASE_DIR"/*/; do
  VCF_FILE="${VCF_DIR}/filtered.PASS.vcf.gz"
  
  # Check if the VCF file exists in this directory
  if [[ -f "$VCF_FILE" ]]; then
    OUT_VCF="${VCF_DIR}/filtered.PASS.chr.vcf.gz"

    echo "Renaming chromosomes for VCF file: $VCF_FILE"
    
    # Rename chromosomes using the mapping file
    bcftools annotate --rename-chrs "$chr_map_file" -Oz -o "$OUT_VCF" "$VCF_FILE"
    
    # Index the renamed VCF file
    tabix -f -p vcf "$OUT_VCF"

    echo "Renaming completed for $VCF_FILE. Output saved to $OUT_VCF"
  else
    echo "No filtered.PASS.vcf.gz found in $VCF_DIR, skipping..."
  fi
done

echo "Chromosome renaming process completed."
