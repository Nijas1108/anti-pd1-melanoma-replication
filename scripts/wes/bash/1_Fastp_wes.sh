#!/bin/bash
#SBATCH --job-name=fastp_trim
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/fastp_%j.out
#SBATCH --error=logs/fastp_%j.err

module purge; module load bluebear
module load bear-apps/2023a
module load fastp/0.24.0-GCC-12.3.0

## decide what samples you want to run this for vanderbilt or ucla
#default vanderbilt

#UCLA Samples
#INPUT_DIR="/rds/projects/e/elhamsak-group5/main_project/data/UCLA_samples" 
#OUT_DIR="/rds/projects/e/elhamsak-group5/main_project/processed_data/trimmed_files/UCLA_samples"

#Vanderbilt Samples
INPUT_DIR="/rds/projects/e/elhamsak-group5/main_project/data/vanderbilt_samples"
OUT_DIR="/rds/projects/e/elhamsak-group5/main_project/processed_data/trimmed_files/vanderbilt_samples"

#UCLA
#SRR_LIST=(
#SRR3083837 SRR3083838 SRR3083839 SRR3083840 SRR3083841 SRR3083842
#SRR3083843 SRR3083844 SRR3083845 SRR3083846 SRR3083847 SRR3083848
#SRR3083849 SRR3083850 SRR3083851 SRR3083852 SRR3083853 SRR3083854
#SRR3083855 SRR3083856 SRR3083857 SRR3083858 SRR3083859 SRR3083860
#SRR3083861 SRR3083862 SRR3083863 SRR3083864 SRR3083866 SRR3083867
#SRR3083868 SRR3083869 SRR3083870 SRR3083871 SRR3083872 SRR3083873
#SRR3083874 SRR3083875 SRR3083876 SRR3083877 SRR3083878 SRR3083879
#SRR3083880 SRR3083881 SRR3083882 SRR3083883 SRR4289714 SRR4289715
#SRR4289716 SRR4289717 SRR4289718 SRR4289719 SRR4289720 SRR4289721
#SRR4289722 SRR4289723 SRR4289724 SRR4289725 SRR4289726 )

##Vanderbilt
SRR_LIST=(
SRR4289728 SRR4289730 SRR4289732 SRR4289734 SRR4289736 SRR4289738
SRR4289740 SRR4289744 SRR4289727 SRR4289729 SRR4289731 SRR4289733
SRR4289735 SRR4289737 SRR4289739 SRR4289743)

for SRR in "${SRR_LIST[@]}"; do
  echo "=============================="
  echo "Checking if trimmed files exist for $SRR"
  echo "=============================="

  OUT_R1="${OUT_DIR}/${SRR}_1.trimmed.fastq.gz"
  OUT_R2="${OUT_DIR}/${SRR}_2.trimmed.fastq.gz"

  # Check if both trimmed files exist
  if [[ -f "$OUT_R1" && -f "$OUT_R2" ]]; then
    echo "Trimmed files already exist for $SRR. Skipping..."
  else
    echo "Running fastp for $SRR"

    R1="${INPUT_DIR}/${SRR}_1.fastq.gz"
    R2="${INPUT_DIR}/${SRR}_2.fastq.gz"

    HTML="${OUT_DIR}/${SRR}.fastp.html"
    JSON="${OUT_DIR}/${SRR}.fastp.json"

    fastp \
      -i "$R1" -I "$R2" \
      -o "$OUT_R1" -O "$OUT_R2" \
      -h "$HTML" -j "$JSON" \
      -w 8 # uses 8 threads

    echo "Done: $SRR"
  fi
done

echo "All fastp jobs completed."
