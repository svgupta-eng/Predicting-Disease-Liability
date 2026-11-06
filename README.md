# Predicting Disease Liability: Chromosome 22 eQTL Mapping & Polygenic Risk Score (PRS) Pipeline

This repository contains the computational pipeline for single-gene and chromosome-wide *cis*-eQTL mapping, p-value threshold tuning across deterministic sample splits (80% Train / 10% Tune / 10% Test), evaluating out-of-sample gene expression prediction (R^2), and calculating Height Polygenic Risk Scores (PRS) using 1000 Genomes reference genotypes.

---

## 1. Environment & Software Dependencies

### System Binaries
* **PLINK 2** (`v2.0.0-b.1-dev`): The PLINK 2 binary is located inside `./plink2_folder/plink2`. Ensure execution permissions are granted before running:
  ```bash
  chmod +x ./plink2_folder/plink2

  Python Dependencies
Python Version: 3.10+

Required Libraries: pandas, numpy, scipy, statsmodels, matplotlib, seaborn

Setup Environment
Bash
pip install pandas numpy scipy statsmodels matplotlib seaborn jupyterlab

2. Data Access & Directory Layout
All reference datasets, genotype files, and expression matrices are placed directly in the repository root directory.

File Sources
1000 Genomes European Reference Panel (Chr 22): 1000G.EUR.22.bed, 1000G.EUR.22.bim, 1000G.EUR.22.fam, 1000G.EUR.1_freqs.afreq, 1000G.EUR.22_freq.afreq, 1000G.EUR.22.GeneX_raw.raw, 1000G.EUR.22.GeneX.bim

GeUVADIS RNA-Seq Expression Data: GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz

Gene Annotations: gene_annot.txt.gz

Height GWAS Summary Statistics: Height.gwas.txt.gz

Repository Layout

Predicting-Disease-Liability/
├── plink2_folder/
│   └── plink2
├── 1000G.EUR.1_freqs.afreq
├── 1000G.EUR.22_freq.afreq
├── 1000G.EUR.22.bed
├── 1000G.EUR.22.bim
├── 1000G.EUR.22.fam
├── 1000G.EUR.22.GeneX.bim
├── 1000G.EUR.22.GeneX_raw.raw
├── GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz
├── Height.gwas.txt.gz
├── gene_annot.txt.gz
├── chr22_protein_coding_eQTLs.txt
├── train_summary_stats.txt
├── train_80.txt
├── tune_10.txt
├── test_10.txt
├── checkpoint.ipynb
└── README.md

3. Reproduction & Execution Instructions
Follow these steps to reproduce all data processing, threshold tuning, and evaluation:

Step 1: Clone Repository
Bash
git clone [https://github.com/svgupta-eng/Predicting-Disease-Liability.git](https://github.com/svgupta-eng/Predicting-Disease-Liability.git)
cd Predicting-Disease-Liability
Step 2: Grant Binary Permissions
Bash
chmod +x ./plink2_folder/plink2
./plink2_folder/plink2 --version
Step 3: Run the Pipeline
Launch JupyterLab and execute the analysis notebook (checkpoint.ipynb) sequentially from top to bottom:

Bash
jupyter lab
Key Workflow Stages:
Sample Alignment & Splitting: Intersects overlapping sample IDs between 1000 Genomes and GeUVADIS data, creating deterministic train_80.txt (80%), tune_10.txt (10%), and test_10.txt (10%) splits.

cis-eQTL Mapping: Fits Ordinary Least Squares (OLS) models for variants within a +/- 500 kb window around target coordinates.

Threshold Tuning: Calculates polygenic scores across p-value thresholds (1, 0.05, 0.01, 0.001, 1e-05) using PLINK 2 with --read-freq to select the threshold maximizing tuning R^2.

Out-of-Sample Test Evaluation: Evaluates out-of-sample R^2 on held-out test samples with variance guardrails for monomorphic test variants.

Chromosome 22 Genome-Wide Loop: Loops over all protein-coding genes on Chromosome 22 and writes results to chr22_protein_coding_eQTLs.txt.

GWAS Polygenic Risk Scoring: Filters Height GWAS summary statistics (p <= 0.0001) and constructs population PRS distributions.

4. Key Output Files
chr22_protein_coding_eQTLs.txt: Tab-separated summary table containing Target Gene, Top SNP, minimum p-value, optimal threshold, and test set R^2.

train_summary_stats.txt: OLS effect sizes and p-values calculated on the 80% training set.

prs_tune_*.sscore & prs_test_final.sscore: Individual polygenic scores output by PLINK 2 for tuning and final testing steps.
