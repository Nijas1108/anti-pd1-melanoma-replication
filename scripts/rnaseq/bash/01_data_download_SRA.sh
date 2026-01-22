#!/bin/bash
#SBATCH --job-name=rnaseq_download
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=rnaseq_%j.out
#SBATCH --error=rnaseq_%j.err

# load SRA Toolkit module on BlueBEAR
module purge
module load bear-apps/2023a
module load SRA-Toolkit/3.0.10-gompi-2023a

# go to your project directory
cd /rds/projects/e/elhamsak-genomics-ngs/group_5/main_project/data/rnaseq

# list your SRA accessions here
SRR_LIST=(
"SRR3184279"
"SRR3184280"
"SRR3184281"
"SRR3184282"
"SRR3184283"
"SRR3184284"
"SRR3184285"
"SRR3184286"
"SRR3184287"
"SRR3184288"
"SRR3184289"
"SRR3184290"
"SRR3184291"
"SRR3184292"
"SRR3184293"
"SRR3184294"
"SRR3184295"
"SRR3184296"
"SRR3184297"
"SRR3184298"
"SRR3184299"
"SRR3184300"
"SRR3184301"
"SRR3184302"
"SRR3184303"
"SRR3184304"
"SRR3184305"
"SRR3184306"
)

for SRR in "${SRR_LIST[@]}"; do
    echo "Processing $SRR ..."

    prefetch $SRR

    fasterq-dump $SRR --split-files --threads 8

    gzip ${SRR}_1.fastq
    gzip ${SRR}_2.fastq

    echo "$SRR done!"
    echo "------------------------"
done

echo "All SRR accessions processed."