# Genomic and Transcriptomic Analysis of Response to Anti–PD-1 Therapy

This repository contains a class-based replication and extension of the study:

**“Genomic and Transcriptomic Features of Response to Anti–PD-1 Therapy in Metastatic Melanoma.”**

The project reproduces and explores key RNA-seq and whole-exome sequencing (WES) analyses used to investigate molecular correlates of response to immune checkpoint blockade therapy. All analyses were performed in a high-performance computing (HPC) environment using reproducible bioinformatics workflows.

---

## Project Objectives

- Reproduce key transcriptomic and genomic analyses from a landmark melanoma immunotherapy study
- Implement an end-to-end RNA-seq and WES analysis pipeline
- Explore molecular features associated with response to anti–PD-1 therapy
- Practice reproducible research and HPC-based bioinformatics workflows
- Create a reusable framework for future extension and research

---

## Repository Structure

```text
.
├── scripts/
│   ├── rnaseq/            # RNA-seq processing & DE analysis
│   │   ├── bash/          # FASTQ → counts pipeline
│   │   └── R/             # DESeq2 & downstream analysis
│   └── wes/
│       ├── bash/          # Alignment, variant calling, annotation
│       └── R/             # WES downstream analyses
│
├── figures/
│   ├── rnaseq/            # RNA-seq plots
│   └── wes/               # WES plots
│
├── presentation/
│   └── class_presentation.pdf
│
├── .gitignore
└── README.md
```

---

## Analysis Overview

### RNA-seq Workflow
1. Data download from public repositories
2. Quality control using **fastp**
3. Alignment using **STAR**
4. Gene quantification via **featureCounts**
5. Differential expression analysis using **DESeq2**
6. Visualization and downstream interpretation

### WES Workflow
1. Read preprocessing and QC
2. Alignment to hg38
3. Variant calling using **GATK Mutect2**
4. Variant filtering and annotation using **Funcotator**
5. Downstream mutation analysis and summarization

---

## Results

This repository contains **summary-level outputs and visualizations only**.

Due to size and data governance considerations, the following are **not included**:
- Raw FASTQ files
- BAM/VCF files
- Reference genomes
- Intermediate alignment files

All results can be reproduced using the provided scripts.

For an overview of the analysis, methodology, and interpretation of results, please refer to the presentation included in this repository:

`presentation/class_presentation.pdf`

This presentation summarizes the biological context, analytical approach, and key findings derived from the RNA-seq and WES analyses.

---

## Data Availability

Raw data were obtained from publicly available datasets referenced in the original publication.

Due to size constraints, raw sequencing files are **not included** in this repository.

Download scripts are provided in:
scripts/rnaseq/bash/
scripts/wes/bash/

---

## Reproducibility

This project was developed and executed on an HPC environment.

### Key tools used:
- fastp
- STAR
- featureCounts
- GATK (Mutect2, Funcotator)
- R (DESeq2, ggplot2)
- CIBERSORTx 
- Bash scripting

The repository is structured to allow full reproduction of results given access to the raw data.

---

## Contributors

This project was completed as a collaborative class assignment.

- **Aishwarya Padmanaban** — WES pipeline development and variant interpretation  
- **Asta Perl** — WES pipeline development and variant interpretation  
- **Aman Kumar** — RNA-seq pipeline development and analysis

All contributors participated equally in the design and execution of this project.

---

## Notes

- This repository was created as part of an MSc-level genomics and bioinformatics course.
- The project may be extended in future work to include:
  - Integrated multi-omics analysis
  - Machine learning–based response prediction
  - Additional melanoma cohorts

---

## Citation

If using or referencing this work, please cite the original study:

> Hugo et al., Cell, 2016. Genomic and Transcriptomic Features of Response to Anti–PD-1 Therapy in Metastatic Melanoma
> https://www.cell.com/cell/fulltext/S0092-8674(16)30215-X?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS009286741630215X%3Fshowall%3Dtrue

---

## Disclaimer

This repository is intended for educational and research purposes only.  
No clinical decisions should be made based on the contents of this repository.
