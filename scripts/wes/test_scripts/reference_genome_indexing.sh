#!/bin/bash
#SBATCH --job-name=ref_genome_indexing
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/bwa_%j.out
#SBATCH --error=logs/bwa_%j.err

module purge; module load bluebear
module load bear-apps/2023a
module load BWA/0.7.17-GCCcore-12.3.0

cd /rds/projects/e/elhamsak-group5/main_project/reference_genome
bwa index hg38.fa
