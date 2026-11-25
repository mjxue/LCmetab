# Code for: High-dimensional serum metabolomic profiling identifies prognostic signatures and host-environment metabolic associations in lung cancer

Maojie Xue, Yuening Ding, Jiacheng Zhou, et al.

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
 
## 2. Installation guide

1. Ensure **R version 4.2.0** is installed on your system.
2. Open the R console and install the required packages using the following command:
   ```R
   install.packages(c("tidyverse", "data.table", "gtsummary", "foreach", "doParallel", 
                      "ggplot2", "ggalluvial", "ggpattern", "patchwork", "circlize", 
                      "ggrepel", "ggbeeswarm", "scales", "survival", "rms", "bnlearn", 
                      "CMAverse", "timeROC", "metamisc", "survIDINRI", "readxl", "plyr"))
