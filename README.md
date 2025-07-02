# BRCA_github

This repository contains scripts and analysis workflows for processing and analyzing BRCA whole-genome sequencing (WGS) data — from raw FASTQ files to somatic variant calling and downstream analysis.

---

## 📁 Repository Structure

```

BRCA_github/
├── scripts/
│   ├── \[scripts to run somatic tools: Brass, Gridss, Manta, Caveman...]
│   └── create\_scripts/
│       └── \[scripts for QC, alignment, variant calling, etc.]
├── data\_analysis/
│   └── \[downstream analysis of somatic variants, visualization, annotation, etc.]

````

---

## 🔬 Workflow Overview

1. **Quality Control** (FastQC, TrimGalore)
2. **Alignment** (BWA-MEM to reference genome)
3. **Post-processing** (MarkDuplicates, BQSR)
4. **Variant Calling**
   - Germline: GATK HaplotypeCaller
   - Somatic: Caveman, Strelka, Manta, Brass, Gridss
5. **Filtering and Annotation**
6. **Data Analysis**
   - Coverage checks
   - Variant frequency and distribution
   - Custom filtering and statistical summaries

---
```

### 2. Execution

Each script `.sh` will create a `.sbatch` file to be run on slurm. For example:

```bash
chmod +x scripts/create_scripts/01_QualityControl.sh
sbatch destination_path/01_QualityControl.sbatch
```

---

## 🛠 Requirements

* SLURM
* BWA, GATK, SAMtools
* Manta, Strelka, Caveman, Brass, Gridss
* TrimGalore, FastQC, Qualimap

---

## 📊 Data Analysis

Post-calling analyses, including variant filtering, summary statistics, and visualization, are in `data_analysis/`.

