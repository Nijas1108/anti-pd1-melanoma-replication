#!/bin/bash
#SBATCH --job-name=STAR_index
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=star_index_%j.out
#SBATCH --error=star_index_%j.err

cd /rds/projects/e/elhamsak-group5/main_project/reference_genome/STAR_index

module purge
module load bear-apps/2024a
module load STAR/2.7.11b-GCC-13.3.0

STAR \
  --runMode genomeGenerate \
  --runThreadN 8 \
  --genomeDir /rds/projects/e/elhamsak-group5/main_project/reference_genome/STAR_index \
  --genomeFastaFiles /rds/projects/e/elhamsak-group5/main_project/reference_genome/hg38.fa \
  --sjdbGTFfile /rds/projects/e/elhamsak-group5/main_project/reference_genome/STAR_index/gencode.v49.annotation.nochr.gtf \
  --sjdbOverhang 99