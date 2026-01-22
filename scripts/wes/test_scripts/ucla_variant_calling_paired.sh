#!/bin/bash
#SBATCH --job-name=mutect2_ucla_pairs
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --array=0-1
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

module purge
module load bluebear
module load bear-apps/2022b
module load GATK/4.4.0.0-GCCcore-12.2.0-Java-17
module load SAMtools/1.17-GCC-12.2.0
module load BCFtools/1.17-GCC-12.2.0

mkdir -p logs

# -------------------------------
# Inputs
# -------------------------------
REF=/rds/projects/e/elhamsak-group5/main_project/reference_genome/hg38.fa
GNOMAD=/rds/projects/e/elhamsak-group5/main_project/reference_genome/gnomAD/af-only-gnomad.hg38.nochr.vcf.gz
COMMON=/rds/projects/e/elhamsak-group5/main_project/reference_genome/gnomAD/small_exac_common_3.hg38.nochr.vcf.gz

ALIGN_DIR=/rds/projects/e/elhamsak-group5/main_project/aligned_files/UCLA_samples
OUT_BASE=/rds/projects/e/elhamsak-group5/main_project/variant_calling/UCLA_pairs/rerun11-01-26

# Optional but strongly recommended for WES speed.
# If you have a BED for the capture kit, put it here.
# If you DON'T have one yet, set BED="" and the script will run genome-wide (slow).
BED="" 
# Example if you do have one:
# BED=/rds/projects/e/elhamsak-group5/main_project/reference_genome/Twist_Exome_Core_hg38_padded.bed

# -------------------------------
# Pairs 
# -------------------------------
# TUMORS=(
# SRR3083866 SRR3083839 SRR3083841 SRR3083837 SRR3083857
# SRR3083863 SRR3083870 SRR4289715 SRR3083849 SRR4289717
# SRR4289719 SRR3083882 SRR3083855 SRR3083868 SRR3083845
# SRR4289721 SRR4289723 SRR3083859 SRR3083847 SRR3083872
# SRR3083861 SRR3083878 SRR3083874 SRR3083853 SRR3083880
# SRR3083843 SRR4289725 SRR4289726 SRR3083851 SRR3083876
# )

# NORMALS=(
# SRR3083867 SRR3083840 SRR3083842 SRR3083838 SRR3083858
# SRR3083864 SRR3083871 SRR4289714 SRR3083850 SRR4289716
# SRR4289718 SRR3083883 SRR3083856 SRR3083869 SRR3083846
# SRR4289720 SRR4289722 SRR3083860 SRR3083848 SRR3083873
# SRR3083862 SRR3083879 SRR3083875 SRR3083854 SRR3083881
# SRR3083844 SRR4289724 SRR4289724 SRR3083852 SRR3083877
# )

TUMORS=(
SRR4289725 SRR4289726
)
# SRR3083866 SRR4289721

NORMALS=(
SRR4289724 SRR4289724
)
# SRR3083867 SRR4289720

TUMOR=${TUMORS[$SLURM_ARRAY_TASK_ID]}
NORMAL=${NORMALS[$SLURM_ARRAY_TASK_ID]}

TUMOR_BAM=${ALIGN_DIR}/${TUMOR}.sorted.bam
NORMAL_BAM=${ALIGN_DIR}/${NORMAL}.sorted.bam

OUTDIR=${OUT_BASE}/${TUMOR}vs${NORMAL}
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "Pair: ${TUMOR} vs ${NORMAL}"
echo "Tumor BAM: $TUMOR_BAM"
echo "Normal BAM: $NORMAL_BAM"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

# -------------------------------
# Reference indexes
# -------------------------------
[ -f "${REF}.fai" ] || samtools faidx "$REF"
[ -f "${REF%.fa}.dict" ] || gatk CreateSequenceDictionary -R "$REF" -O "${REF%.fa}.dict"

# Threading: helps GKL/PairHMM actually use the CPUs you requested
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_THREAD_LIMIT=${SLURM_CPUS_PER_TASK}

# Intervals for Mutect2/FilterMutectCalls (WES speedup)
INTERVAL_ARGS=()
if [[ -n "${BED}" ]]; then
  [[ -f "${BED}" ]] || { echo "BED file not found: ${BED}"; exit 1; }
  INTERVAL_ARGS=(-L "${BED}")
  echo "Using WES targets BED: ${BED}"
else
  echo "No BED provided: running genome-wide (slow)."
fi

# -------------------------------
# Step 1: Mutect2 (tumor-normal)
# -------------------------------
if [[ ! -f unfiltered.vcf.gz ]]; then
  echo "Running Mutect2..."
  gatk --java-options "-Xmx40G" Mutect2 \
    -R "$REF" \
    -I "$TUMOR_BAM" -tumor "$TUMOR" \
    -I "$NORMAL_BAM" -normal "$NORMAL" \
    --germline-resource "$GNOMAD" \
    --f1r2-tar-gz f1r2.tar.gz \
    --native-pair-hmm-threads "$SLURM_CPUS_PER_TASK" \
    "${INTERVAL_ARGS[@]}" \
    -O unfiltered.vcf.gz
else
  echo "Skipping Mutect2 (unfiltered.vcf.gz exists)."
fi

if [[ ! -f unfiltered.vcf.gz.tbi ]]; then
  gatk IndexFeatureFile -I unfiltered.vcf.gz
fi

# -------------------------------
# Step 2: Orientation bias model
# -------------------------------
if [[ ! -f read-orientation-model.tar.gz ]]; then
  echo "Learning read orientation model..."
  gatk LearnReadOrientationModel \
    -I f1r2.tar.gz \
    -O read-orientation-model.tar.gz
else
  echo "Skipping LearnReadOrientationModel (exists)."
fi

# -------------------------------
# Step 3: Pileups + contamination
# -------------------------------
if [[ ! -f pileups.table ]]; then
  echo "Running GetPileupSummaries..."
  gatk GetPileupSummaries \
    -I "$TUMOR_BAM" \
    -V "$COMMON" \
    -L "$COMMON" \
    -O pileups.table
else
  echo "Skipping GetPileupSummaries (exists)."
fi

if [[ ! -f contamination.table ]]; then
  echo "Running CalculateContamination..."
  gatk CalculateContamination \
    -I pileups.table \
    -O contamination.table \
    --tumor-segmentation segments.table
else
  echo "Skipping CalculateContamination (exists)."
fi

# -------------------------------
# Step 4: Filter calls
# -------------------------------
if [[ ! -f filtered.vcf.gz ]]; then
  echo "Running FilterMutectCalls..."
  gatk FilterMutectCalls \
    -R "$REF" \
    -V unfiltered.vcf.gz \
    --contamination-table contamination.table \
    --tumor-segmentation segments.table \
    --ob-priors read-orientation-model.tar.gz \
    "${INTERVAL_ARGS[@]}" \
    -O filtered.vcf.gz
else
  echo "Skipping FilterMutectCalls (filtered.vcf.gz exists)."
fi

# -------------------------------
# Step 5: Keep PASS only
# -------------------------------
if [[ ! -f filtered.PASS.vcf.gz ]]; then
  echo "Extracting PASS variants..."
  bcftools view -f PASS filtered.vcf.gz -Oz -o filtered.PASS.vcf.gz
  bcftools index -t filtered.PASS.vcf.gz
else
  echo "Skipping PASS extraction (exists)."
fi

echo "Done for ${TUMOR} vs ${NORMAL}"
ls -lh unfiltered.vcf.gz filtered.vcf.gz filtered.PASS.vcf.gz || true