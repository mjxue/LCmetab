rm(list=ls())
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(survival)

setwd('D:/Cancer Model/03LC/LC_Metab/v3/')

tissue_raw <- read_xlsx('01data/04Evolutionary metabolic landscape from preneoplasia to invasive lung adenocarcinoma/Supplementary Data 2.xlsx', sheet = 'tissue sample ')

############################ QC procedure#########################################
###################### 1. Cohort1_tissue##########################################
metablist <- tissue_raw$Metabolite
data_matrix <- tissue_raw[, -1]

data_t <- as.matrix(t(data_matrix))
colnames(data_t) <- metablist

sample_names <- colnames(data_matrix)
sample_info <- data.frame(
  sample = sample_names,
  type = case_when(
    grepl("QC", sample_names) ~ "QC",
    grepl("_N", sample_names) ~ "N",
    grepl("_T", sample_names) ~ "T",
    TRUE ~ "unknown"
  ),
  batch = gsub("(Batch\\d+).*", "\\1", sample_names),
  id = str_extract(sample_names, "(?<=_N|_T|_QC_)\\d+")
)

qc_idx <- sample_info$type == "QC"
test_idx <- sample_info$type %in% c("N", "T")

qc_data <- data_t[qc_idx, , drop = FALSE]
qc_all_mean <- colMeans(qc_data, na.rm = TRUE)

batches <- unique(sample_info$batch[sample_info$type == "QC"])
batch_qc_means <- matrix(NA, nrow = length(batches), ncol = ncol(data_t))
rownames(batch_qc_means) <- batches
colnames(batch_qc_means) <- colnames(data_t)

for (b in batches) {
  idx_b <- which(sample_info$batch == b & sample_info$type == "QC")
  if (length(idx_b) > 0) {
    batch_qc_means[b, ] <- colMeans(data_t[idx_b, , drop = FALSE], na.rm = TRUE)
  }
}

data_corrected <- data_t
for (i in which(test_idx)) {
  batch_i <- sample_info$batch[i]
  qc_adj <- batch_qc_means[batch_i, ]
  factor <- ifelse(qc_adj == 0 | is.na(qc_adj), 1, qc_all_mean / qc_adj)
  data_corrected[i, ] <- data_t[i, ] * factor
}


tumor_idx <- sample_info$type == "T"
data_tumor <- data_corrected[tumor_idx, , drop = FALSE]
# Exclude Batch01_T090 where the signal is abnormally too low.
data_tumor <- data_tumor[rownames(data_tumor) != "Batch01_T090", ]

tumor_row_sums <- rowSums(data_tumor, na.rm = TRUE)
data_tumor_normalized <- sweep(data_tumor, 1, tumor_row_sums, "/")

# log transformation
tumor_info <- sample_info[sample_info$type %in% c("T"), ]

data_log <- as.data.frame(log(data_tumor_normalized))
data_scaled <- as.data.frame(scale(data_log, center = TRUE, scale = TRUE))

data_scaled$sampleID <- substr(rownames(data_scaled), 10, nchar(rownames(data_scaled)))

cohort1_clin <- read_xlsx('01data/04Evolutionary metabolic landscape from preneoplasia to invasive lung adenocarcinoma/Supplementary Data 1.xlsx', sheet = 'Patients in cohort 1')
# cohort1_tissue
cohort1_clin <- cohort1_clin %>%
  mutate(
    Ctype  = 'LUAD',
    Cstage = ifelse(Stage == '0',0,
                    ifelse(Stage %in% c('Ⅰ','I'),1,
                           ifelse(Stage %in% c('Ⅱ','II'),2,
                                  ifelse(Stage %in% c('Ⅲ'),3,4)))),
    Gender = ifelse(Gender=='Male',0,1),
    sampleID = substr(cohort1_clin$`Patient ID`, 2, nchar(cohort1_clin$`Patient ID`)),
    Smoke = `Smoke (0: Non-smoker;1: Smoker)`
  ) %>%
  filter(Cstage != 0)

cohort1_data <- left_join(cohort1_clin, data_scaled)

rmlist <- c('batch_qc_means','data_corrected','data_log','data_matrix','data_tumor_normalized','data_scaled','data_t','tumor_info',
            'b','batch_i','batches','factor','i','idx_b','qc_adj','qc_all_mean','qc_idx','tumor_row_sums','sample_names','sample_info',
            'test_idx','qc_data','tumor_idx')
rm(list = rmlist)


metablist <- tissue_raw$Metabolite

# Auxiliary function: Determines whether backticks are required
needs_backtick <- function(x) {
  grepl("[^a-zA-Z0-9_.]", x)  
}

results <- data.frame(
  metab = character(),
  HR = numeric(),
  CIL = numeric(),
  CIU = numeric(),
  Pval = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_along(metablist)) {
  metab <- metablist[i]
  
  formula_str <- paste0('Surv(OS, OS_status) ~ Age + as.factor(Cstage) + Smoke + Gender + `', metab, '`')
  
  fit <- tryCatch(
    coxph(as.formula(formula_str), data = cohort1_data),
    error = function(e) NULL
  )
  
  if (!is.null(fit)) {
    sum <- summary(fit)
    coef_names <- rownames(sum$coefficients)
    
    if (metab %in% coef_names) {
      coef_name <- metab
    } else if (paste0("`", metab, "`") %in% coef_names) {
      coef_name <- paste0("`", metab, "`")
    } else {
      coef_name <- coef_names[gsub("`", "", coef_names) == metab][1]
    }
    
    if (!is.na(coef_name)) {
      results <- rbind(results, data.frame(
        metab = metab,
        HR = sum$coefficients[coef_name, "exp(coef)"],
        CIL = sum$conf.int[coef_name, "lower .95"],
        CIU = sum$conf.int[coef_name, "upper .95"],
        Pval = sum$coefficients[coef_name, "Pr(>|z|)"],
        stringsAsFactors = FALSE
      ))
    } else {
      results <- rbind(results, data.frame(
        metab = metab,
        HR = NA,
        CIL = NA,
        CIU = NA,
        Pval = NA,
        stringsAsFactors = FALSE
      ))
    }
  } else {
    results <- rbind(results, data.frame(
      metab = metab,
      HR = NA,
      CIL = NA,
      CIU = NA,
      Pval = NA,
      stringsAsFactors = FALSE
    ))
  }
}

head(results)

write.csv(results, '03result/02Table/TabS12_LUAD_valid.csv',row.names = F)

###### TabS9 #################################################
combined_data <- cohort1_data %>%
  select(c('Patient ID','sampleID','OS','OS_status','Age','Cstage','Smoke','Gender')) %>%
  mutate(

    group = case_when(
      Cstage == 1 ~ 'Stage I',
      Cstage == 2 ~ 'Stage II',
      Cstage == 3 ~ 'Stage III',
      Cstage == 4 ~ 'Stage IV',
      TRUE ~ NA_character_
    ),
    Dead = factor(OS_status,
                  levels = c(0,1),
                  labels = c('Alive', 'Dead')),
    Gender = factor(Gender,
                    levels = c(0,1),
                    labels = c('Male', 'Female')),
    Smoke = factor(Smoke,
                   levels = c(0, 1, 2),
                   labels = c("Never", "Former", "Current"))
  ) %>%
  filter(!is.na(group)) 


baseline_variables <- c('Age', "Gender", "Smoke", 'Dead')


table_one <- combined_data %>%
  select(all_of(baseline_variables), group) %>%
  tbl_summary(
    by = group, # 按 4 个分期进行分组
    type = all_dichotomous() ~ "categorical",
    label = list(
      Age ~ "Age at diagnosis (years)",
      Gender ~ "Gender",
      Smoke ~ "Smoke status",
      Dead ~ "Vital status"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})", 
      all_categorical() ~ "{n} ({p})"     
    ),
    digits = list(
      all_continuous() ~ c(2, 2),         
      all_categorical() ~ c(0, 2)         
    ),
    missing_text = "Unknown"
  ) %>%
  add_overall(last = TRUE) %>%
  modify_header(
    label = "Characteristic",
    all_stat_cols() ~ "**{level}**\nN = {n}"
  )

table_one


table_one_df <- table_one %>%
  as_tibble()


write.csv(table_one_df,
          file = "03result/02Table/TabS9.csv",
          row.names = FALSE,
          fileEncoding = "UTF-8") 
