################ Table1 calculate ##############################
rm(list = ls())
library(readxl)
library(dplyr)
library(survival)
library(gtsummary)

setwd('D:/Cancer Model/03LC/01LC_Metab/v3/')
load('01data/02cleandata/metab_compute.RData')
load('01data/02cleandata/data_origin.RData')
# source('02code/function/Biomarker_screening_function.R')

metab_compute$smoke_bin <- ifelse(metab_compute$Smoke==0,0,ifelse(is.na(metab_compute$Smoke),NA,1))

expose_list <- data_origin$Index[data_origin$Origin=='Exogenous']
metab_list <- data_origin$Index[data_origin$Origin=='Endogenous']

covar_list <- c('Age_diagnosis','Age_blood',"Gender","Smoke","Pack_year","Ctype","Cstage",'radio','chemo','smoke_bin','bl2di','timediff')
metab_compute <- subset(metab_compute,select = c("ID",covar_list,'OS_diag','OS_bld','Dead','PFS_diag','PFS',data_origin$Index))

combined_data <- metab_compute %>%
  mutate(
    group = case_when(
      Cstage == 1 ~ 'Stage I',
      Cstage == 2 ~ 'Stage II',
      Cstage == 3 ~ 'Stage III',
      Cstage == 4 ~ 'Stage IV',
      TRUE ~ NA_character_
    ),
    Dead = factor(Dead,
                  levels = c(0,1),
                  labels = c('Alive', 'Dead')),
    Gender = factor(Gender,
                    levels = c(0,1),
                    labels = c('Male', 'Female')),
    Smoke = factor(Smoke,
                   levels = c(0, 1, 2),
                   labels = c("Never", "Former", "Current")),
    Ctype = factor(Ctype,
                   levels = c('LUAD', 'LUSC', 'Other_NSCLC', 'Other'),
                   labels = c("LUAD", "LUSC", "Other NSCLC", "Other")),
    radio = factor(radio,
                   levels = c(0, 1),
                   labels = c("No", "Yes")),
    chemo = factor(chemo,
                   levels = c(0,1),
                   labels = c('No','Yes'))
  ) %>%
  filter(!is.na(group))

baseline_variables <- c('Age_diagnosis', 'Age_blood', "Gender", "Smoke", "Ctype", "radio",'chemo', 'Dead')

table_one <- combined_data %>%
  select(all_of(baseline_variables), group) %>%
  tbl_summary(
    by = group, # 按 4 个分期进行分组
    type = all_dichotomous() ~ "categorical",
    label = list(
      Age_diagnosis ~ "Age at diagnosis (years)",
      Age_blood ~ "Age at blood sampling (years)",
      Gender ~ "Gender (%)",
      Smoke ~ "Smoke status (%)",
      Ctype ~ "Histological subtype (%)",
      radio ~ "Radiotherapy (%)",
      chemo ~ "Chemotherapy (%)",
      Dead ~ "Vital status (%)"
    ),
    statistic = list(
      all_continuous() ~ "{mean} +- {sd}", 
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
          file = "03result/02Table/Tab1.csv",
          row.names = FALSE,
          fileEncoding = "UTF-8")

