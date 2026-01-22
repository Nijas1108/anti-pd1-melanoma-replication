#!/bin/bash
#SBATCH --job-name=gnomad_v3_download
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=logs/gnomad_v3_download_%j.out
#SBATCH --error=logs/gnomad_v3_download_%j.err

# Load modules
module purge
module load bluebear
module load bear-apps/2023a
module load gcloud/472.0.0
module load bear-apps/2024a
module load BCFtools/1.21-GCC-13.3.0

# Directories
OUT_DIR=/rds/projects/e/elhamsak-group5/main_project/reference_genome/gnomAD
mkdir -p "$OUT_DIR"
mkdir -p logs

# Download per-chromosome genome VCFs and indexes
for chr in {1..22} X Y; do
    echo "Downloading chromosome $chr"
    gsutil cp \
        gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr${chr}.vcf.bgz \
        "$OUT_DIR/"
    gsutil cp \
        gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr${chr}.vcf.bgz.tbi \
        "$OUT_DIR/"
done

# Merge all chromosomes into a single VCF
echo "Merging all chromosomes into one VCF"
bcftools concat -a -O z \
    "$OUT_DIR"/gnomad.genomes.v3.1.2.sites.chr*.vcf.bgz \
    -o "$OUT_DIR"/gnomad.genomes.v3.1.2.merged.vcf.gz

# Index the merged VCF
echo "Indexing merged VCF"
bcftools index "$OUT_DIR"/gnomad.genomes.v3.1.2.merged.vcf.gz

echo "Merged gnomAD v3.1.2 genome VCF is ready at $OUT_DIR/gnomad.genomes.v3.1.2.merged.vcf.gz"
