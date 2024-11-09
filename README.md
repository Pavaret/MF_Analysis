# Thai Macaca fascicularis Genomic Analysis

This repository provides scripts for generating figures and conducting analysis presented in the manuscript:

**Genomic Characterization of Thai Macaca fascicularis based on PacBio HiFi Sequencing**

## Table of Contents
- [Introduction](#introduction)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Data Access](#data-access)
- [Running the Analysis](#running-the-analysis)
- [Results Overview](#results-overview)

---

## Introduction

This project offers a genomic insight into the Thai *Macaca fascicularis* using high-fidelity long-read sequencing technology (PacBio HiFi). Comparative analyses with the Cambodian *M. fascicularis* reference genome reveal SNP and SV patterns, which underscore the genomic diversity within this species.

## Repository Structure

The repository includes:
- **Scripts**:
  - `Variant_Analyser.sh`: Main bash script to run the analysis pipeline, calling Python and R scripts in sequence.
  - `main.py`: Python script that organizes functions for generating figures from data files.
  - `Figure1.py`, `Figure2.py`, `Figure3.py`, `Figure4.py`: Python scripts for generating individual figures.
  - `Figure5.R`: R script for functional enrichment analysis in Figure 5.

## Requirements

Install required libraries and dependencies:
- **Python 3.8+**
  - `pandas`
  - `numpy`
  - `matplotlib`
  - `scipy`
  - `seaborn`
  - `altair` (for Vega-Altair visualization)
- **R**
  - `ggplot2`
  - `dplyr`
  - `gprofiler2`

You can install Python dependencies with:
```bash
pip install pandas numpy matplotlib scipy seaborn altair
```

For R, install packages by running:
```R
install.packages(c("ggplot2", "dplyr", "gprofiler2"))
```

## Data Access

The sequencing data is available in the NCBI SRA repository under [PRJNA1077753](https://www.ncbi.nlm.nih.gov/sra/PRJNA1077753). Download this data to reproduce the analyses.

## Running the Analysis

To reproduce the figures described in the manuscript, you can use `Variant_Analyser.sh`, which integrates and runs all the individual scripts in sequence.

```bash
bash Variant_Analyser.sh
```

This command will execute each of the scripts to produce Figures 1 through 5.

### Individual Script Execution
If you prefer to run individual scripts separately, execute the following commands:

```bash
# Figure 1: Genetic variant distribution
python Figure1.py

# Figure 2: SNP density across chromosomes
python Figure2.py

# Figure 3: Structural variant types
python Figure3.py

# Figure 4: SV distribution across chromosomes
python Figure4.py

# Figure 5: Functional enrichment analysis
Rscript Figure5.R
```

## Results Overview

The analyses reveal critical insights into the genomic architecture of Thai *M. fascicularis*, including:
- SNP and SV distributions across chromosomal regions.
- Enrichment of structural variants affecting cellular, neuronal, and reproductive processes.
- Comparative patterns highlighting genetic uniqueness in Thai populations of *M. fascicularis*.

---
