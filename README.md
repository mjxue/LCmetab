# Code for: High-dimensional serum metabolomic profiling identifies prognostic signatures and host-environment metabolic associations in lung cancer

Maojie Xue, Yuening Ding, Jiacheng Zhou, et al.

> **Data Availability Notice**  
> Due to strict patient privacy and data use agreements, the raw data from the Boston Lung Cancer Study (BLCS) and the UK Biobank (UKB) cannot be publicly shared. However, to fulfill the code availability and reproducibility requirements, a **de-identified subset dataset (`demo_BLCS_data.RData`)** is provided in this repository for testing the core scripts. The provided R scripts are fully functional and are designed to accept standard data formats (e.g., `.RData` or `.csv`) containing survival outcomes, clinical covariates, and metabolic features.

This repository contains the R source code used for the statistical analyses, biomarker screening, sensitivity analyses, multi-omics integration, and figure generation in the manuscript.

## 1. System requirements
### Software dependencies
- **R version 4.2.0** (The code was developed and tested on this version).
- **Required R Packages**: 
  The following R packages and versions were used in our environment. 
  - `tidyverse` (2.0.0)
  - `data.table` (1.16.2)
  - `gtsummary` (2.3.0)
  - `foreach` (1.5.2)
  - `doParallel` (1.0.17)
  - `ggplot2` (4.0.3)
  - `ggalluvial` (0.12.5)
  - `ggpattern` (1.3.1)
  - `patchwork` (1.3.0)
  - `circlize` (0.4.16)
  - `ggrepel` (0.9.6)
  - `ggbeeswarm` (0.7.2)
  - `scales` (1.4.0)
  - `survival` (3.6.4)
  - `rms` (6.8.0)
  - `bnlearn` (4.8.1)
  - `CMAverse` (0.1.0)
  - `timeROC` (0.4)
  - `metamisc` (0.4.0)
  - `survIDINRI` (1.1.2)
  - `readxl` (1.4.3)
  - `plyr` (1.8.9)
  - `survminer` (0.5.0)

### Hardware requirements
- No non-standard hardware is required. The software can be run on a standard desktop computer.

## 2. Installation guide

1. Ensure **R version 4.2.0** is installed on your system.
2. Open the R console and install the required packages using the following command:
   ```R
   install.packages(c("tidyverse", "data.table", "gtsummary", "foreach", "doParallel", 
                      "ggplot2", "ggalluvial", "ggpattern", "patchwork", "circlize", 
                      "ggrepel", "ggbeeswarm", "scales", "survival", "rms", "bnlearn", 
                      "CMAverse", "timeROC", "metamisc", "survIDINRI", "readxl", "plyr",
                      "survminer"))
   ```
   
3.Typical install time: Installation should take approximately 10-20 minutes on a normal desktop computer, depending on your internet connection and compilation time.

## 3. Demo
As stated in the Data Availability notice, the raw data from the BLCS and UKB cannot be publicly shared due to patient privacy restrictions. To fulfill the code availability requirement and allow reviewers to verify our pipeline, we have provided a **de-identified subset dataset with positive results (`demo_BLCS_data.RData`)**. This dataset specifically supports the execution of `13Fig4_plot.R`, which evaluates the predictive performance (AUC, E/O ratio, DCA, and NRI/IDI) of the metabolite-augmented model versus the baseline clinical model.

### Instructions to run on data:
1. **Download required files:** Ensure the following files are in your working directory:
   - `demo_BLCS_data.RData` (contains the demo_BLCS_data data frame).
   - `stdca.R` (the custom function required for Decision Curve Analysis, which is available in the `function/` folder of this repository).
2. **Modify the paths:** Open the script `13Fig4_plot.R`. 
   - **CRITICAL:** Comment out or delete the `setwd('D:/Cancer Model/...')` command.
   - **CRITICAL:** Modify the line `source('D:/Cancer Model/my_function/stdca.R')` to point to your local path where `stdca.R` is located (e.g., `source('Function/stdca.R')`).
3. **Create output directories:** Ensure the folders `03result/01Figure/` and `03result/02Table/` exist in your working directory. If not, create them before running the script (e.g., `dir.create("03result/01Figure", recursive = TRUE)`).
4. **Run the script:** Execute the entire `13Fig4_plot.R` script. It will automatically build the Cox models, perform 10-fold cross-validation, calculate the metrics, and generate the plots/tables.

### Expected output:
Upon successful execution, the script will generate and save the following files in your working directory:
- **Figures (saved in `03result/01Figure/`):**
  - `Fig4A_AUC.pdf`: Bar chart showing the time-dependent AUC of the two models from 1 to 10 years.
  - `Fig4B_EO.pdf`: Plot showing the Expected/Observed (E/O) ratio for model calibration.
  - `Fig4C_DCA.pdf`: Four-panel Decision Curve Analysis (DCA) plots for 3-year and 5-year survival.
- **Tables (saved in `03result/02Table/`):**
  - `TabS13_AUC_compare_optbia.csv`: Table comparing AUCs with optimism-corrected values.
  - `EO_compare.csv`: Table of E/O ratios at different time points.
  - `TabS14_IDI_NRI_1to10yr.csv`: Table containing the calculated NRI and IDI metrics for 1 to 10 years.
- **Console output:** The R console will print progress messages (e.g., "Processing time = X years") during the IDI/NRI calculation loop, and print the final `all_results` table.

### Expected run time for demo:
- On a "normal" desktop computer (e.g., 4-core CPU, 16GB RAM), running `13Fig4_plot.R` with the demo dataset will take approximately **5-15 minutes**. The majority of this time is consumed by the 1000-iteration bootstrap in the `IDI.INF` function and the 10-fold cross-validation.

## 4. Instructions for use
### File Structure
- 01Creat_Tab1.R to 11Sub_analysis_plot.R, and 13Fig4_plot.R: Scripts for the BLCS metabolomic analysis and visualization.
- 12External_validation.R: Scripts for the Chinese lung cancer cohort analyses.
- 14UKBPPP_clinical.R and 15multi_Fig5.R: Scripts for the UKB-PPP trans-omics analyses.

- Function/: Contains custom R functions (e.g., Biomarker_screening_function.R, stdca.R) essential for the analysis. Ensure paths are correctly sourced in the main scripts.

## 5. How to run on your data
1. Clone or download this repository to your local machine.

2. Modify the file paths in the scripts (specifically the setwd() and load() commands) to point to your local directories and datasets.

3. Adjust the n_cores parameter in scripts utilizing doParallel (e.g., 02BIF_permutation.R and 03Biomarker_screen.R) to match your computer's available CPU cores.
