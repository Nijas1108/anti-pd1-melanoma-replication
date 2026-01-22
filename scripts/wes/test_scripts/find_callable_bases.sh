#!/bin/bash
#SBATCH --job-name=mosdepth_tmb
#SBATCH --account elhamsak-group5
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --array=0-38
#SBATCH --output=logs/mosdepth_%A_%a.out
#SBATCH --error=logs/mosdepth_%A_%a.err

# 1. Load Environment
module purge
module load bear-apps/2022b
module load Miniforge3/24.1.2-0
module load BEDTools/2.30.0-GCC-12.2.0
source activate mosdepth_env

# 2. Define Lists
TUMOUR_LIST=(
SRR3083866 SRR3083839 SRR3083841 SRR3083837 SRR3083857 SRR3083863 SRR3083870 SRR4289715 SRR3083849 SRR4289717
SRR4289719 SRR3083882 SRR3083855 SRR3083868 SRR3083845 SRR4289721 SRR4289723 SRR3083859 SRR3083847 SRR3083872
SRR3083861 SRR3083878 SRR3083874 SRR3083853 SRR3083880 SRR3083843 SRR4289725 SRR4289726 SRR3083851 SRR3083876
SRR4289728 SRR4289730 SRR4289732 SRR4289734 SRR4289736 SRR4289738 SRR4289740 SRR4289742 SRR4289744
)

NORMAL_LIST=(
SRR3083867 SRR3083840 SRR3083842 SRR3083838 SRR3083858 SRR3083864 SRR3083871 SRR4289714 SRR3083850 SRR4289716
SRR4289718 SRR3083883 SRR3083856 SRR3083869 SRR3083846 SRR4289720 SRR4289722 SRR3083860 SRR3083848 SRR3083873
SRR3083862 SRR3083879 SRR3083875 SRR3083854 SRR3083881 SRR3083844 SRR4289724 SRR4289724 SRR3083852 SRR3083877
SRR4289727 SRR4289729 SRR4289731 SRR4289733 SRR4289735 SRR4289737 SRR4289739 SRR4289741 SRR4289743
)

# 3. Clean and Extract IDs
# Using xargs to trim any potential invisible whitespace
T_ID=$(echo "${TUMOUR_LIST[$SLURM_ARRAY_TASK_ID]}" | xargs)
N_ID=$(echo "${NORMAL_LIST[$SLURM_ARRAY_TASK_ID]}" | xargs)

TARGET_BED="/rds/projects/e/elhamsak-group5/main_project/reference_genome/Twist_Exome_Core_hg38_padded.bed"
OUT_DIR="/rds/projects/e/elhamsak-group5/main_project/aligned_files/coverage_results"
mkdir -p "$OUT_DIR" logs

# 4. Find BAM Paths
UCLA_PATH="/rds/projects/e/elhamsak-group5/main_project/aligned_files/UCLA_samples/processed_bams"
VAND_PATH="/rds/projects/e/elhamsak-group5/main_project/aligned_files/vanderbilt_samples/processed_bams"

if [[ -f "${UCLA_PATH}/${T_ID}.final.bam" ]]; then
    T_BAM="${UCLA_PATH}/${T_ID}.final.bam"
    N_BAM="${UCLA_PATH}/${N_ID}.final.bam"
elif [[ -f "${VAND_PATH}/${T_ID}.final.bam" ]]; then
    T_BAM="${VAND_PATH}/${T_ID}.final.bam"
    N_BAM="${VAND_PATH}/${N_ID}.final.bam"
else
    echo "ERROR: Files for index $SLURM_ARRAY_TASK_ID ($T_ID / $N_ID) not found."
    exit 1
fi

echo "Processing Index $SLURM_ARRAY_TASK_ID: Tumor $T_ID vs Normal $N_ID"

# check to see if mosdepth env is loaded properly
which mosdepth
mosdepth --version

# 5. Run Mosdepth (Using single quotes for quantize string)
# mosdepth -t 4 -n --by "$TARGET_BED" --quantize '0:10:' "${OUT_DIR}/${T_ID}" "$T_BAM"
# mosdepth -t 4 -n --by "$TARGET_BED" --quantize '0:10:' "${OUT_DIR}/${N_ID}" "$N_BAM"

# Correct quantize usage
mosdepth -t 4 --by "$TARGET_BED" "${OUT_DIR}/${T_ID}" "$T_BAM"
mosdepth -t 4 --by "$TARGET_BED" "${OUT_DIR}/${N_ID}" "$N_BAM"


# 6. Extract 10x regions and Intersect
# We use -f to force zcat to work if file exists
# zcat "${OUT_DIR}/${T_ID}.quantize.bed.gz" | awk '$4 >= 10' > "${OUT_DIR}/${T_ID}_T_10x.bed"
# zcat "${OUT_DIR}/${N_ID}.quantize.bed.gz" | awk '$4 >= 10' > "${OUT_DIR}/${N_ID}_N_10x.bed"

zcat "${OUT_DIR}/${T_ID}.regions.bed.gz" | awk '$4 >= 10' > "${OUT_DIR}/${T_ID}_T_10x.bed"
zcat "${OUT_DIR}/${N_ID}.regions.bed.gz" | awk '$4 >= 10' > "${OUT_DIR}/${N_ID}_N_10x.bed"


bedtools intersect -a "${OUT_DIR}/${T_ID}_T_10x.bed" -b "${OUT_DIR}/${N_ID}_N_10x.bed" > "${OUT_DIR}/${T_ID}_joint.bed"

# 7. Calculate Megabases
TOTAL_BASES=$(awk '{sum += ($3 - $2)} END {print (sum == "" ? 0 : sum)}' "${OUT_DIR}/${T_ID}_joint.bed")
TOTAL_MB=$(echo "scale=4; $TOTAL_BASES / 1000000" | bc)

echo -e "${T_ID}\t${N_ID}\t${TOTAL_BASES}\t${TOTAL_MB}" > "${OUT_DIR}/${T_ID}_summary.txt"

# 8. Clean up
# rm "${OUT_DIR}/${T_ID}_T_10x.bed" "${OUT_DIR}/${N_ID}_N_10x.bed" "${OUT_DIR}/${T_ID}_joint.bed"

# run after the job is complete: 
# cd /rds/projects/e/elhamsak-group5/main_project/aligned_files/coverage_results
# cat *_summary.txt > all_samples_callable_bases.txt
