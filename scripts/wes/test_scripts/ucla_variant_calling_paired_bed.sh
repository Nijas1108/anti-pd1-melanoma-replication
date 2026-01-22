#!/bin/bash
#SBATCH --job-name=mutect2_paired
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=0-1
#SBATCH --output=logs/mutect2_%A_%a.out
#SBATCH --error=logs/mutect2_%A_%a.err

set -euo pipefail

# --------------------
# Load modules
# --------------------
module purge
module load bluebear
module load bear-apps/2022b
module load GATK/4.4.0.0-GCCcore-12.2.0-Java-17
module load SAMtools/1.17-GCC-12.2.0
module load BCFtools/1.17-GCC-12.2.0

# --------------------
# Paths
# --------------------
REF_FASTA=/rds/projects/e/elhamsak-group5/main_project/reference_genome/hg38.fa
BAM_DIR=/rds/projects/e/elhamsak-group5/main_project/aligned_files/UCLA_samples
BED_FILE=/rds/projects/e/elhamsak-group5/main_project/reference_genome/Twist_Exome_Core_hg38_padded.bed
GNOMAD_VCF=/rds/projects/e/elhamsak-group5/main_project/reference_genome/gnomAD/af-only-gnomad.hg38.nochr.vcf.gz
OUT_DIR=/rds/projects/e/elhamsak-group5/main_project/variants/UCLA_samples/paired_bed_extra_filters/

mkdir -p logs "$OUT_DIR"

# --------------------
# Sample lists (ORDER MATTERS)
# --------------------
TUMOUR_LIST=(
SRR3083866 SRR3083839 SRR3083841 SRR3083837 SRR3083857
SRR3083863 SRR3083870 SRR4289715 SRR3083849 SRR4289717
SRR4289719 SRR3083882 SRR3083855 SRR3083868 SRR3083845
SRR4289721 SRR4289723 SRR3083859 SRR3083847 SRR3083872
SRR3083861 SRR3083878 SRR3083874 SRR3083853 SRR3083880
SRR3083843 SRR4289725 SRR4289726 SRR3083851 SRR3083876
)

NORMAL_LIST=(
SRR3083867 SRR3083840 SRR3083842 SRR3083838 SRR3083858
SRR3083864 SRR3083871 SRR4289714 SRR3083850 SRR4289716
SRR4289718 SRR3083883 SRR3083856 SRR3083869 SRR3083846
SRR4289720 SRR4289722 SRR3083860 SRR3083848 SRR3083873
SRR3083862 SRR3083879 SRR3083875 SRR3083854 SRR3083881
SRR3083844 SRR4289724 SRR4289724 SRR3083852 SRR3083877
)

# --------------------
# Select pair
# --------------------
TUMOUR=${TUMOUR_LIST[$SLURM_ARRAY_TASK_ID]}
NORMAL=${NORMAL_LIST[$SLURM_ARRAY_TASK_ID]}

TUMOUR_BAM=${BAM_DIR}/${TUMOUR}.sorted.bam
NORMAL_BAM=${BAM_DIR}/${NORMAL}.sorted.bam

[[ -f "$TUMOUR_BAM" ]] || { echo "Missing tumour BAM"; exit 1; }
[[ -f "$NORMAL_BAM" ]] || { echo "Missing normal BAM"; exit 1; }


PREFIX=${OUT_DIR}/${TUMOUR}_vs_${NORMAL}

echo "Tumour: $TUMOUR"
echo "Normal: $NORMAL"

# --------------------
# 1) Mutect2 (unfiltered + F1R2)
# --------------------
if [[ ! -f ${PREFIX}.unfiltered.vcf.gz ]]; then
    echo "Running Mutect2..."
    gatk --java-options "-Xmx48G" Mutect2 \
        -R "$REF_FASTA" \
        -I "$TUMOUR_BAM" -tumor "$TUMOUR" \
        -I "$NORMAL_BAM" -normal "$NORMAL" \
        -L "$BED_FILE" \
        --native-pair-hmm-threads 8 \
        --f1r2-tar-gz "${PREFIX}.f1r2.tar.gz" \
        -O "${PREFIX}.unfiltered.vcf.gz"
fi

# --------------------
# 2) Pileup summaries
# --------------------
gatk GetPileupSummaries \
    -I "$NORMAL_BAM" \
    -V "$GNOMAD_VCF" \
    -L "$BED_FILE" \
    -O "${PREFIX}.normal.pileup.table"

gatk GetPileupSummaries \
    -I "$TUMOUR_BAM" \
    -V "$GNOMAD_VCF" \
    -L "$BED_FILE" \
    -O "${PREFIX}.tumour.pileup.table"

# --------------------
# 3) Contamination
# --------------------
gatk CalculateContamination \
    -I "${PREFIX}.tumour.pileup.table" \
    -matched "${PREFIX}.normal.pileup.table" \
    -O "${PREFIX}.contamination.table"

# --------------------
# 4) Read orientation model
# --------------------
gatk LearnReadOrientationModel \
    -I "${PREFIX}.f1r2.tar.gz" \
    -O "${PREFIX}.read-orientation-model.tar.gz"

# --------------------
# 5) Filter calls
# --------------------
gatk FilterMutectCalls \
    -R "$REF_FASTA" \
    -V "${PREFIX}.unfiltered.vcf.gz" \
    --contamination-table "${PREFIX}.contamination.table" \
    --ob-priors "${PREFIX}.read-orientation-model.tar.gz" \
    -L "$BED_FILE" \
    -O "${PREFIX}.filtered.vcf.gz"

# --------------------
# 6) Index filtered VCF
# --------------------
gatk IndexFeatureFile \
    -I "${PREFIX}.filtered.vcf.gz"

# --------------------
# 7) Extract PASS
# --------------------
bcftools view \
    -f PASS \
    "${PREFIX}.filtered.vcf.gz" \
    -Oz -o "${PREFIX}.filtered.PASS.vcf.gz"

bcftools index "${PREFIX}.filtered.PASS.vcf.gz"

echo "Done: ${PREFIX}.filtered.PASS.vcf.gz"