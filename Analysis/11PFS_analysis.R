rm(list = ls())

# 加载必要包
library(dplyr)
library(survival)

# 设置工作目录（请根据实际路径调整）
setwd('D:/Cancer Model/03LC/LC_Metab/v3/')

load('01data/02cleandata/metab_compute.RData')   
load('01data/02cleandata/data_origin.RData')    
load('03result/Ultra_Robust_Biomarkers.RData')   


source('02code/function/Biomarker_screening_function.R')

# main effect
main_indices <- c(robust_exp_main$Index, robust_metab_main$Index)
# age-inter
inter_indices <- c(robust_exp_age$Index, robust_metab_age$Index)

#  define covar ---------------------------------------------
covars_diag <- c("Gender", "Smoke", "Ctype", "Cstage", "Age_blood")

#  outcome PFS COX --------------------------------------------

# main effect
pfs_main <- full_cox_main(
  data        = metab_compute,
  feature_list = main_indices,
  covars      = covars_diag,
  time_col    = "PFS_diag",
  event_col   = "PFS"
)

# age-inter
pfs_inter <- full_cox_inter(
  data           = metab_compute,
  feature_list   = inter_indices,
  interaction_var = "Age_blood",
  covars         = covars_diag,
  time_col       = "PFS_diag",
  event_col      = "PFS"
)


anno <- data_origin %>% select(Index, Compounds, Origin)

pfs_main_anno <- pfs_main %>%
  left_join(anno, by = "Index") %>%
  select(Index, Compounds, Origin, HR, LCI, UCI, P_Value) %>%
  rename(HR_PFS=HR,LCI_PRS=LCI, UCI_PFS=UCI, P_Value_PFS=P_Value) %>%
  arrange(P_Value_PFS) %>%
  left_join(dplyr::bind_rows(
    robust_exp_main,
    robust_metab_main)
    )

pfs_inter_anno <- pfs_inter %>%
  left_join(anno, by = "Index") %>%
  select(Index, Compounds, Origin, HR, LCI, UCI, P_Value) %>%
  rename(HR_PFS=HR,LCI_PRS=LCI, UCI_PFS=UCI, P_Value_PFS=P_Value) %>%
  arrange(P_Value_PFS) %>%
  left_join(dplyr::bind_rows(
    robust_exp_age,
    robust_metab_age)
  )


write.csv(pfs_main_anno, file = "03result/02Table/TabS10.csv", row.names = FALSE)
write.csv(pfs_inter_anno, file = "03result/02Table/TabS11.csv", row.names = FALSE)

