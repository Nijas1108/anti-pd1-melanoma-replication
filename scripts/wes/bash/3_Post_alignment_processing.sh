#!/bin/bash
#SBATCH --job-name=preproc
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=0-15
#SBATCH --output=logs/prep_%a.out

set -euo pipefail

## decide what samples you want to run this for vanderbilt or ucla
#default vanderbilt
#uncomment paths for UCLA and the sample list if you want to use it

#Load modules
module purge
module load bluebear
module load bear-apps/2022b
module load GATK/4.4.0.0-GCCcore-12.2.0-Java-17
module load SAMtools/1.17-GCC-12.2.0

# All Paths
#General
REF="/rds/projects/e/elhamsak-group5/main_project/reference_genome/hg38/hg38.fa"
KNOWN_SITES="/rds/projects/e/elhamsak-group5/main_project/reference_genome/gnomAD/"

# UCLA
#IN_DIR="/rds/projects/e/elhamsak-group5/main_project/processed_data/aligned_files/UCLA_samples"
#OUT_DIR="/rds/projects/e/elhamsak-group5/main_project/processed_data/aligned_files/UCLA_samples/processed_bams"

# Vanderbilt
IN_DIR="/rds/projects/e/elhamsak-group5/main_project/processed_data/aligned_files/vanderbilt_samples"
OUT_DIR="/rds/projects/e/elhamsak-group5/main_project/processed_data/aligned_files/vanderbilt_samples/processed_bams"

mkdir -p "$OUT_DIR"

#Sample Lists
#UCLA 
#SAMPLES=(
#SRR3083837 SRR3083838 SRR3083839 SRR3083840 SRR3083841
#SRR3083842 SRR3083843 SRR3083844 SRR3083845 SRR3083846
#SRR3083847 SRR3083848 SRR3083849 SRR3083850 SRR3083851 
#SRR3083852 SRR3083853 SRR3083854 SRR3083855 SRR3083856
#SRR3083857 SRR3083858 SRR3083859 SRR3083860 SRR3083861
#SRR3083862 SRR3083863 SRR3083864 SRR3083866 SRR3083867
#SRR3083868 SRR3083869 SRR3083870 SRR3083871 SRR3083872
#SRR3083873 SRR3083874 SRR3083875 SRR3083876 SRR3083877
#SRR3083878 SRR3083879 SRR3083880 SRR3083881 SRR3083882
#SRR3083883 SRR4289714 SRR4289715 SRR4289716 SRR4289717
#SRR4289718 SRR4289719 SRR4289720 SRR4289721 SRR4289722
#SRR4289723 SRR4289724 SRR4289725 SRR4289726
#)

#Vanderbilt
SAMPLES=(
SRR4289728 SRR4289730 SRR4289732 SRR4289734 SRR4289736 SRR4289738
SRR4289740 SRR4289744 SRR4289727 SRR4289729 SRR4289731 SRR4289733
SRR4289735 SRR4289737 SRR4289739 SRR4289743
)

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

echo "Processing Sample: $SAMPLE"

# ---------------------------------------------------------
# STEP 1: Mark Duplicates 
# ---------------------------------------------------------
# REMOVE_DUPLICATES=false (marks them so Mutect2 can ignore them)
gatk --java-options "-Xmx32G" MarkDuplicates \
    -I "${IN_DIR}/${SAMPLE}.sorted.bam" \
    -O "${OUT_DIR}/${SAMPLE}.marked.bam" \
    -M "${OUT_DIR}/${SAMPLE}_dup_metrics.txt" \
    --CREATE_INDEX true \
    --VALIDATION_STRINGENCY SILENT
# --validation stringency- 3 levels, silent- no warnings. 

# Base Quality Score Recalibration (BSQR)
# ---------------------------------------------------------
# STEP 2: BQSR Phase 1 (Model Building)
# ---------------------------------------------------------
gatk --java-options "-Xmx32G" BaseRecalibrator \
    -I "${OUT_DIR}/${SAMPLE}.marked.bam" \
    -R "$REF" \
    --known-sites "$KNOWN_SITES" \
    -O "${OUT_DIR}/${SAMPLE}_recal_data.table"

# ---------------------------------------------------------
# STEP 3: BQSR Phase 2 (Applying Correction)
# ---------------------------------------------------------
gatk --java-options "-Xmx32G" ApplyBQSR \
    -I "${OUT_DIR}/${SAMPLE}.marked.bam" \
    -R "$REF" \
    --bqsr-recal-file "${OUT_DIR}/${SAMPLE}_recal_data.table" \
    -O "${OUT_DIR}/${SAMPLE}.final.bam"

# ---------------------------------------------------------
# STEP 4: Cleanup
# ---------------------------------------------------------
# Keep the .final.bam and .final.bai; delete the intermediate .marked.bam
rm "${OUT_DIR}/${SAMPLE}.marked.bam"
rm "${OUT_DIR}/${SAMPLE}.marked.bai"

echo "Workflow complete for $SAMPLE. Ready for Variant Calling with Mutect2."
