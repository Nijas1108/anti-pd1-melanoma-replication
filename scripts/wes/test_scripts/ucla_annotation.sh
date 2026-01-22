#!/bin/bash
#SBATCH --job-name=ucla_annotate_vcf
#SBATCH --cpus-per-task=2 # Funcotator is single thread-based?
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=0-29
#SBATCH --output=logs/annotate_vcf_%A_%a.out
#SBATCH --error=logs/annotate_vcf_%A_%a.err

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

echo "Job started at: $(date)"

# --------------------
# Paths
# --------------------
REF_FASTA=/rds/projects/e/elhamsak-group5/main_project/reference_genome/hg38_chr/hg38_chr.fa
IN_DIR=/rds/projects/e/elhamsak-group5/main_project/variant_calling/ucla_processed_bed
OUT_DIR=/rds/projects/e/elhamsak-group5/main_project/annotation/ucla_processed_bed

mkdir -p "$OUT_DIR"

# Funcotator data sources
FUNCOTATOR_DIR=/rds/projects/e/elhamsak-group5/main_project/reference_genome/funcotator/funcotator_dataSources.v1.7.20200521g

# --------------------
# Sample lists (ORDER MATTERS!)
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

PAIR_DIR="${IN_DIR}/${TUMOUR}vs${NORMAL}"
PASS_VCF="${PAIR_DIR}/filtered.PASS.chr.vcf.gz"

PREFIX_OUT=${OUT_DIR}/${TUMOUR}_vs_${NORMAL}

echo "Tumour: $TUMOUR"
echo "Normal: $NORMAL"
echo "Pair dir: $PAIR_DIR"

# --------------------
# Check if PASS VCF exists
# --------------------
if [[ ! -f $PASS_VCF ]]; then
    echo "PASS VCF does not exist. Skipping annotation."
    exit 1
else
    echo "Found PASS VCF: $PASS_VCF"
fi

# --------------------
# 1) Run Funcotator Annotation
# --------------------
if [[ ! -f ${PREFIX_OUT}.annotated.vcf.gz ]]; then
    echo "Running Funcotator..."
    gatk Funcotator \
        -R "$REF_FASTA" \
        -V "$PASS_VCF" \
        --ref-version hg38 \
        --output "${PREFIX_OUT}.annotated.vcf.gz" \
        --data-sources-path "$FUNCOTATOR_DIR" \
        --output-file-format VCF \
        --java-options "-Xmx48G"
else
    echo "Annotated VCF already exists: ${PREFIX_OUT}.annotated.vcf.gz"
fi

# --------------------
# 2) Index annotated VCF
# --------------------
if [[ ! -f "${PREFIX_OUT}.annotated.vcf.gz.tbi" ]]; then
    echo "Indexing annotated VCF..."
    gatk IndexFeatureFile -I "${PREFIX_OUT}.annotated.vcf.gz"
else
    echo "Annotated VCF already indexed."
fi

echo "Annotation completed: ${PREFIX_OUT}.annotated.vcf.gz"

# --------------------
# 3) Validation
# --------------------
gatk ValidateVariants \
  -R "$REF_FASTA" \
  -V ${PREFIX_OUT}.annotated.vcf.gz

echo "Job finished at: $(date)"
