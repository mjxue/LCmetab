################ Table calculate ##############################
rm(list = ls())
library(readxl)
library(dplyr)
library(survival)
library(gtsummary)

setwd('D:/Cancer Model/03LC/LC_Metab/v3/')
load('01data/03protein/ukb_protein_analysis.RData')

combined_data <- ukb_protein_analysis %>%
  mutate(
    Dead = factor(Dead,
                  levels = c(0,1),
                  labels = c('Alive', 'Dead')),
    Gender = factor(sex,
                    levels = c(1,2),
                    labels = c('Male', 'Female')),
    Smoke = factor(smokeingstatus,
                   levels = c(0, 1, 2),
                   labels = c("Never", "Former", "Current")),
    Ctype = factor(histo_new,
                   levels = c('LUAD','LUSC','Other NSCLC', 'SCLC', 'Others'),
                   labels = c('LUAD','LUSC','Other NSCLC', 'SCLC', 'Others'))
  )


baseline_variables <- c('age', "Gender", "Smoke", "Ctype", 'Dead')


table_one <- combined_data %>%
  select(all_of(baseline_variables)) %>%
  tbl_summary(
    type = all_dichotomous() ~ "categorical",
    label = list(
      age ~ "Age at diagnosis (years)",
      Gender ~ "Gender",
      Smoke ~ "Smoke status",
      Ctype ~ "Histological subtype",
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
  # add_overall(last = TRUE) %>%
  modify_header(
    label = "Characteristic",
    all_stat_cols() ~ "**{level}**\nN = {n}"
  )

table_one


table_one_df <- table_one %>%
  as_tibble()

write.csv(table_one_df,
          file = "03result/TabS15_UKBPPP.csv",
          row.names = FALSE,
          fileEncoding = "UTF-8") 

