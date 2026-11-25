##########################################Analysis######################################
rm(list = ls())
library(readxl)
library(dplyr)
library(survival)

setwd('D:/Cancer Model/03LC/01LC_Metab/v3/')
load('01data/02cleandata/metab_compute.RData')
load('01data/02cleandata/data_origin.RData')
source('02code/function/Biomarker_screening_function.R')

metab_compute$smoke_bin <- ifelse(metab_compute$Smoke==0,0,ifelse(is.na(metab_compute$Smoke),NA,1))

expose_list <- data_origin$Index[data_origin$Origin=='Exogenous']
metab_list <- data_origin$Index[data_origin$Origin=='Endogenous']

covar_list <- c('Age_diagnosis','Age_blood',"Gender","Smoke","Pack_year","Ctype","Cstage",'radio','chemo','smoke_bin','bl2di','timediff')
metab_compute <- subset(metab_compute,select = c("ID",covar_list,'OS_diag','OS_bld','Dead','PFS_diag','PFS',data_origin$Index))

################### Analysis A: Bootstrap Robustness Analysis ##################################
library(dplyr)
library(foreach)
library(doParallel)
library(survival)

n_cores <- 25
cl <- makeCluster(n_cores)
registerDoParallel(cl)

n_iterations <- 1000 
covars_diag <- c("Gender", "Smoke", "Ctype", "Cstage", "Age_blood")

system.time({  
  results <- foreach(
    i = 1:n_iterations, 
    .packages = c("survival", "dplyr", "caret"), 
    .export = c("metab_compute", "expose_list", "metab_list", "covars_diag", 
                "maincox_analysis", "GEinter_analysis") 
  ) %dopar% {
    
    set.seed(3100 + i)     
    
    # --- A. Bootstrap sampling based on the entire dataset ---
    n_total <- nrow(metab_compute)
    
    boot_index <- sample(1:n_total, size = n_total, replace = TRUE)
    
    train_data <- metab_compute[boot_index, ]
    
    oob_index <- setdiff(1:n_total, unique(boot_index))
    valid_data <- metab_compute[oob_index, ]    
    # --- B. Extract Subsets ---
    expose_train <- subset(train_data, select = c(expose_list, covars_diag, 'OS_diag', 'Dead'))
    expose_valid <- subset(valid_data, select = c(expose_list, covars_diag, 'OS_diag', 'Dead'))
    
    metab_train <- subset(train_data, select = c(metab_list, covars_diag, 'OS_diag', 'Dead'))
    metab_valid <- subset(valid_data, select = c(metab_list, covars_diag, 'OS_diag', 'Dead'))
    
    # --- C1. main effect ---
    diag_exp_main <- maincox_analysis(expose_train, expose_valid, expose_list, covars_diag, time_col = "OS_diag")
    diag_metab_main <- maincox_analysis(metab_train, metab_valid, metab_list, covars_diag, time_col = "OS_diag")
    
    # --- D1. select robust biomarker with main effect ---
    valid_diag_exp_main <- diag_exp_main[
      diag_exp_main$fdr_train <= 0.05 & 
        diag_exp_main$P_Value_valid <= 0.05 & 
        ((diag_exp_main$HR_train > 1 & diag_exp_main$HR_valid > 1) | 
           (diag_exp_main$HR_train < 1 & diag_exp_main$HR_valid < 1)), 
    ]
    
    valid_diag_metab_main <- diag_metab_main[
      diag_metab_main$fdr_train <= 0.05 & 
        diag_metab_main$P_Value_valid <= 0.05 & 
        ((diag_metab_main$HR_train > 1 & diag_metab_main$HR_valid > 1) | 
           (diag_metab_main$HR_train < 1 & diag_metab_main$HR_valid < 1)), 
    ]
    
    # --- C2. matabolite-age interaction ---
    
    diag_exp_age <- GEinter_analysis(expose_train, expose_valid, expose_list, 
                                     interaction_var = "Age_blood", covars = covars_diag)
    diag_metab_age <- GEinter_analysis(metab_train, metab_valid, metab_list, 
                                       interaction_var = "Age_blood", covars = covars_diag)
    
    # --- D2. select robust metabolite with metabolite-age interaction ---
    valid_diag_exp_age <- diag_exp_age[
      diag_exp_age$fdr_train <= 0.05 & 
        diag_exp_age$P_Value_valid <= 0.05 & 
        ((diag_exp_age$HR_train > 1 & diag_exp_age$HR_valid > 1) | 
           (diag_exp_age$HR_train < 1 & diag_exp_age$HR_valid < 1)), 
    ]
    
    valid_diag_metab_age <- diag_metab_age[
      diag_metab_age$fdr_train <= 0.05 & 
        diag_metab_age$P_Value_valid <= 0.05 & 
        ((diag_metab_age$HR_train > 1 & diag_metab_age$HR_valid > 1) | 
           (diag_metab_age$HR_train < 1 & diag_metab_age$HR_valid < 1)), 
    ]
    
    # --- E. Return the 4 types of results from each loop as a list ---
    list(
      exp_main_features   = valid_diag_exp_main$Index,
      metab_main_features = valid_diag_metab_main$Index,
      exp_age_features    = valid_diag_exp_age$Index,
      metab_age_features  = valid_diag_metab_age$Index
    )
  }
})

stopCluster(cl)


robustness_exp_main   <- lapply(results, function(x) x$exp_main_features)
robustness_metab_main <- lapply(results, function(x) x$metab_main_features)
robustness_exp_age    <- lapply(results, function(x) x$exp_age_features)
robustness_metab_age  <- lapply(results, function(x) x$metab_age_features)

# Define an auxiliary function to count frequencies and sort them
get_freq_table <- function(feature_list, col_name) {
  df <- as.data.frame(table(unlist(feature_list)))
  if(nrow(df) > 0) {
    colnames(df) <- c(col_name, "Frequency")
    
    df$Percentage <- paste0(round((df$Frequency / n_iterations) * 100, 1), "%")
    df <- df %>% arrange(desc(Frequency))
  } else {
    df <- data.frame(Feature = character(), Frequency = integer(), Percentage = character())
    colnames(df)[1] <- col_name
  }
  return(df)
}

freq_exp_main   <- get_freq_table(robustness_exp_main, "Exogenous_Main")
freq_metab_main <- get_freq_table(robustness_metab_main, "Endogenous_Main")
freq_exp_age    <- get_freq_table(robustness_exp_age, "Exogenous_Age_Interaction")
freq_metab_age  <- get_freq_table(robustness_metab_age, "Endogenous_Age_Interaction")

save(freq_exp_main, freq_metab_main, freq_exp_age, freq_metab_age, file='03result/freq_summary3100.RData')

################### Analysis B：sensitivity analysis for radiotherapy and chemotherapy ###################
rm(list = ls())
library(readxl)
library(dplyr)
library(survival)

setwd('D:/Cancer Model/03LC/01LC_Metab/v3/')
load('01data/02cleandata/metab_compute.RData')
load('01data/02cleandata/data_origin.RData')
source('02code/function/Biomarker_screening_function.R')

metab_compute$smoke_bin <- ifelse(metab_compute$Smoke==0,0,ifelse(is.na(metab_compute$Smoke),NA,1))

expose_list <- data_origin$Index[data_origin$Origin=='Exogenous']
metab_list <- data_origin$Index[data_origin$Origin=='Endogenous']

covar_list <- c('Age_diagnosis','Age_blood',"Gender","Smoke","Pack_year","Ctype","Cstage",'radio','chemo','smoke_bin','bl2di','timediff')
metab_compute <- subset(metab_compute,select = c("ID",covar_list,'OS_diag','OS_bld','Dead','PFS_diag','PFS',data_origin$Index))

covars_diag <- c("Gender", "Smoke", "Ctype", "Cstage", "Age_blood",'chemo','radio')

# 1. Main
full_exp_main_treat <- full_cox_main(metab_compute, expose_list, covars_diag, time_col = "OS_diag")
full_metab_main_treat <- full_cox_main(metab_compute, metab_list, covars_diag, time_col = "OS_diag")

# 2. age-inter
full_exp_age_treat <- full_cox_inter(metab_compute, expose_list, "Age_blood", covars_diag, time_col = "OS_diag")
full_metab_age_treat <- full_cox_inter(metab_compute, metab_list, "Age_blood", covars_diag, time_col = "OS_diag")

full_exp_main_treat <- full_exp_main_treat %>% dplyr::rename(HR_treat = HR, LCI_treat = LCI, UCI_treat = UCI, P_Value_treat = P_Value)
full_metab_main_treat <- full_metab_main_treat %>% dplyr::rename(HR_treat = HR, LCI_treat = LCI, UCI_treat = UCI, P_Value_treat = P_Value)
full_exp_age_treat <- full_exp_age_treat %>% dplyr::rename(HR_treat = HR, LCI_treat = LCI, UCI_treat = UCI, P_Value_treat = P_Value)
full_metab_age_treat <- full_metab_age_treat %>% dplyr::rename(HR_treat = HR, LCI_treat = LCI, UCI_treat = UCI, P_Value_treat = P_Value)

save(full_exp_main_treat, full_metab_main_treat, full_exp_age_treat, full_metab_age_treat, 
     file='03result/full_alldata_res_treat.RData')

################### Analysis C：delayed cox analysis ###################
rm(list = ls())
library(readxl)
library(dplyr)
library(survival)

setwd('D:/Cancer Model/03LC/01LC_Metab/v3/')
load('01data/02cleandata/metab_compute.RData')
load('01data/02cleandata/data_origin.RData')
source('02code/function/Biomarker_screening_function.R')

metab_compute$smoke_bin <- ifelse(metab_compute$Smoke==0,0,ifelse(is.na(metab_compute$Smoke),NA,1))

expose_list <- data_origin$Index[data_origin$Origin=='Exogenous']
metab_list <- data_origin$Index[data_origin$Origin=='Endogenous']

covar_list <- c('Age_diagnosis','Age_blood',"Gender","Smoke","Pack_year","Ctype","Cstage",'radio','chemo','bl2di','timediff')
metab_compute <- subset(metab_compute,select = c("ID",covar_list,'OS_diag','OS_bld','Dead','PFS_diag','PFS',data_origin$Index))

# ----------  delayed-entry  ----------
metab_compute <- metab_compute %>%
  mutate(
    # 采血早于确诊的，进入时间记为 0；晚于确诊的，记为 bl2di
    entry_time = ifelse(bl2di < 0, 0, bl2di),
    exit_time = OS_diag
  ) %>%
  filter(exit_time > entry_time)

# ---------- covar ----------
covars_diag <- c("Gender", "Smoke", "Ctype", "Cstage", "Age_blood",'chemo','radio')

covars_diag_inter <- setdiff(covars_diag, "Age_blood")


base_covars_main  <- paste(covars_diag,        collapse = " + ")
base_covars_inter <- paste(covars_diag_inter,  collapse = " + ")

# fit Delayed-entry Cox model
# ---------- All features：Exogenous / Endogenous × Main / Interaction ----------
features_to_test <- rbind(
  data.frame(Index = expose_list, Type = "Main",        Origin = "Exogenous", stringsAsFactors = FALSE),
  data.frame(Index = metab_list,  Type = "Main",        Origin = "Endogenous", stringsAsFactors = FALSE),
  data.frame(Index = expose_list, Type = "Interaction", Origin = "Exogenous", stringsAsFactors = FALSE),
  data.frame(Index = metab_list,  Type = "Interaction", Origin = "Endogenous", stringsAsFactors = FALSE)
)

delayed_cox_results <- data.frame()

for (i in seq_len(nrow(features_to_test))) {
  feat      <- features_to_test$Index[i]
  feat_type <- features_to_test$Type[i]
  origin    <- features_to_test$Origin[i]
  
  if (i %% 50 == 0) cat(sprintf("  ... Progress %d / %d (Current: %s, %s)\n",
                                i, nrow(features_to_test), feat, feat_type))
  
  res_row <- tryCatch({
    if (feat_type == "Main") {
      form <- as.formula(paste0(
        "Surv(entry_time, exit_time, Dead) ~ `", feat, "` + ", base_covars_main
      ))
      fit <- coxph(form, data = metab_compute)
      sm  <- summary(fit)
      
      if (!(feat %in% rownames(sm$coefficients))) return(NULL)
      
      data.frame(
        Index = feat, Type = feat_type, Origin = origin,
        HR      = sm$conf.int[feat, "exp(coef)"],
        LCI     = sm$conf.int[feat, "lower .95"],
        UCI     = sm$conf.int[feat, "upper .95"],
        P_Value = sm$coefficients[feat, "Pr(>|z|)"],
        stringsAsFactors = FALSE
      )
    } else {
      form <- as.formula(paste0(
        "Surv(entry_time, exit_time, Dead) ~ `", feat, "` * Age_blood + ",
        base_covars_inter
      ))
      fit <- coxph(form, data = metab_compute)
      sm  <- summary(fit)
      
      term_names <- rownames(sm$coefficients)
      inter_term <- term_names[
        grepl(":", term_names, fixed = TRUE) &
          grepl(feat, term_names, fixed = TRUE)
      ]
      if (length(inter_term) == 0) return(NULL)
      inter_name <- inter_term[1]
      
      data.frame(
        Index = feat, Type = feat_type, Origin = origin,
        HR      = sm$conf.int[inter_name, "exp(coef)"],
        LCI     = sm$conf.int[inter_name, "lower .95"],
        UCI     = sm$conf.int[inter_name, "upper .95"],
        P_Value = sm$coefficients[inter_name, "Pr(>|z|)"],
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) {
    message(sprintf("  [Skip] %s (%s/%s): %s",
                    feat, origin, feat_type, conditionMessage(e)))
    NULL
  })
  
  if (!is.null(res_row)) {
    delayed_cox_results <- rbind(delayed_cox_results, res_row)
  }
}

print(delayed_cox_results)

delayed_exp_main <- delayed_cox_results[delayed_cox_results$Type == "Main" & delayed_cox_results$Origin == 'Exogenous', ]
delayed_metab_main <- delayed_cox_results[delayed_cox_results$Type == "Main" & delayed_cox_results$Origin == 'Endogenous', ]

delayed_exp_age <- delayed_cox_results[delayed_cox_results$Type == "Interaction" & delayed_cox_results$Origin == 'Exogenous', ]
delayed_metab_age <- delayed_cox_results[delayed_cox_results$Type == "Interaction" & delayed_cox_results$Origin == 'Endogenous', ]


delayed_exp_main <- delayed_exp_main %>% dplyr::rename(HR_delay = HR, LCI_delay = LCI, UCI_delay = UCI, P_Value_delay = P_Value)
delayed_metab_main <- delayed_metab_main %>% dplyr::rename(HR_delay = HR, LCI_delay = LCI, UCI_delay = UCI, P_Value_delay = P_Value)
delayed_exp_age <- delayed_exp_age %>% dplyr::rename(HR_delay = HR, LCI_delay = LCI, UCI_delay = UCI, P_Value_delay = P_Value)
delayed_metab_age <- delayed_metab_age %>% dplyr::rename(HR_delay = HR, LCI_delay = LCI, UCI_delay = UCI, P_Value_delay = P_Value)

save(delayed_cox_results, 
     delayed_exp_main, delayed_metab_main, 
     delayed_exp_age, delayed_metab_age,  
     file = "03result/delayed_results.RData")

################################# OS_diag_notreat ############################################

covars_diag <- c("Gender", "Smoke", "Ctype", "Cstage", "Age_blood")
covars_bld <- c("Gender", "Smoke", "Ctype", "Cstage", "Age_blood",'timediff')

# 1. main
full_diag_exp_main <- full_cox_main(metab_compute, expose_list, covars_diag, time_col = "OS_diag")
full_diag_metab_main <- full_cox_main(metab_compute, metab_list, covars_diag, time_col = "OS_diag")

# 2. age-inter
full_diag_exp_age <- full_cox_inter(metab_compute, expose_list, "Age_blood", covars_diag, time_col = "OS_diag")
full_diag_metab_age <- full_cox_inter(metab_compute, metab_list, "Age_blood", covars_diag, time_col = "OS_diag")

full_diag_exp_main <- full_diag_exp_main %>% dplyr::rename(HR_diag = HR, LCI_diag = LCI, UCI_diag = UCI, P_Value_diag = P_Value)
full_diag_metab_main <- full_diag_metab_main %>% dplyr::rename(HR_diag = HR, LCI_diag = LCI, UCI_diag = UCI, P_Value_diag = P_Value)
full_diag_exp_age <- full_diag_exp_age %>% dplyr::rename(HR_diag = HR, LCI_diag = LCI, UCI_diag = UCI, P_Value_diag = P_Value)
full_diag_metab_age <- full_diag_metab_age %>% dplyr::rename(HR_diag = HR, LCI_diag = LCI, UCI_diag = UCI, P_Value_diag = P_Value)

save(full_diag_exp_main, full_diag_metab_main, full_diag_exp_age, full_diag_metab_age, 
     file='03result/full_alldata_diag_res_notreat.RData')


######################################### robust biomarker check ########################################################

rm(list = ls())
setwd('D:/Cancer Model/03LC/01LC_Metab/v3/')
load('01data/02cleandata/data_origin.RData')
load('03result/freq_summary3100.RData')

load('03result/delayed_results.RData')
load('03result/full_alldata_res_treat.RData')
load('03result/full_alldata_diag_res_notreat.RData')
library(dplyr)

# ==========================================
# 0. Set screening criteria
# ==========================================
# Sensitivity analysis P-value threshold
alpha_threshold <- 0.05
# Threshold of Bootstrap Incorporating Frequency (BIF) based on permutation
boot_threshold <- 73.2

# ==========================================
# 1. Endogenous
# ==========================================
# A. Main Effect
robust_metab_main <- freq_metab_main %>%
  rename(Index = Endogenous_Main) %>%
  left_join(full_diag_metab_main %>% select(Index, HR_diag, LCI_diag, UCI_diag, P_Value_diag), by = 'Index') %>%
  left_join(full_metab_main_treat %>% select(Index, HR_treat, LCI_treat, UCI_treat, P_Value_treat), by = "Index") %>%
  left_join(delayed_metab_main %>% select(Index, HR_delay, LCI_delay, UCI_delay, P_Value_delay), by = "Index") %>%
  mutate(Boot_Percent_Num = as.numeric(gsub("%", "", Percentage)))

robust_metab_main <- robust_metab_main %>% filter(
  Boot_Percent_Num >= boot_threshold, 
  P_Value_diag  <= 0.05,
  P_Value_treat <= 0.05,  
  P_Value_delay <= 0.05,        
  (HR_diag>1&HR_treat>1&HR_delay>1)|(HR_diag<1&HR_treat<1&HR_delay<1),
  ) %>%
  arrange(desc(Boot_Percent_Num))


# B. GxE Interaction
robust_metab_age <- freq_metab_age %>%
  rename(Index = Endogenous_Age_Interaction) %>%
  left_join(full_diag_metab_age %>% select(Index, HR_diag, LCI_diag, UCI_diag, P_Value_diag), by = 'Index') %>%
  left_join(full_metab_age_treat %>% select(Index, HR_treat, LCI_treat, UCI_treat, P_Value_treat), by = "Index") %>%
  left_join(delayed_metab_age %>% select(Index, HR_delay, LCI_delay, UCI_delay, P_Value_delay), by = "Index") %>%
  mutate(Boot_Percent_Num = as.numeric(gsub("%", "", Percentage)))

robust_metab_age <- robust_metab_age %>% filter(
  Boot_Percent_Num >= boot_threshold,        
  P_Value_diag  <= 0.05,
  P_Value_treat <= 0.05,                    
  P_Value_delay <= 0.05,                   
  (HR_diag>1&HR_treat>1&HR_delay>1)|(HR_diag<1&HR_treat<1&HR_delay<1),
) %>%
  arrange(desc(Boot_Percent_Num))

# ==========================================
# 2.Exogenous
# ==========================================
# A. Main effect
robust_exp_main <- freq_exp_main %>%
  rename(Index = Exogenous_Main) %>%
  left_join(full_diag_exp_main %>% select(Index, HR_diag, LCI_diag, UCI_diag, P_Value_diag), by = 'Index') %>%
  left_join(full_exp_main_treat %>% select(Index, HR_treat, LCI_treat, UCI_treat, P_Value_treat), by = "Index") %>%
  left_join(delayed_exp_main %>% select(Index, HR_delay, LCI_delay, UCI_delay, P_Value_delay), by = "Index") %>%
  mutate(Boot_Percent_Num = as.numeric(gsub("%", "", Percentage)))

robust_exp_main <- robust_exp_main %>% filter(
  Boot_Percent_Num >= boot_threshold,         
  P_Value_diag  <= 0.05,
  P_Value_treat <= 0.05,                      
  P_Value_delay <= 0.05,                     
  (HR_diag>1&HR_treat>1&HR_delay>1)|(HR_diag<1&HR_treat<1&HR_delay<1),       
) %>%
  arrange(desc(Boot_Percent_Num)) # 按稳健性从高到低排序

# B. 效
robust_exp_age <- freq_exp_age %>%
  rename(Index = Exogenous_Age_Interaction) %>%
  left_join(full_diag_exp_age %>% select(Index, HR_diag, LCI_diag, UCI_diag, P_Value_diag), by = 'Index') %>%
  left_join(full_exp_age_treat %>% select(Index, HR_treat, LCI_treat, UCI_treat, P_Value_treat), by = "Index") %>%
  left_join(delayed_exp_age %>% select(Index, HR_delay, LCI_delay, UCI_delay, P_Value_delay), by = "Index") %>%
  mutate(Boot_Percent_Num = as.numeric(gsub("%", "", Percentage)))

robust_exp_age <- robust_exp_age %>% filter(
  Boot_Percent_Num >= boot_threshold,        
  P_Value_diag  <= 0.05,
  P_Value_treat <= 0.05,                      
  P_Value_delay <= 0.05,                     
  (HR_diag>1&HR_treat>1&HR_delay>1)|(HR_diag<1&HR_treat<1&HR_delay<1),         
) %>%
  arrange(desc(Boot_Percent_Num)) 

save(robust_metab_main, robust_metab_age, robust_exp_main, robust_exp_age, 
     file = "03result/Ultra_Robust_Biomarkers.RData")

################### TabS1#############################################
freq_main <- dplyr::bind_rows(
  freq_exp_main %>% rename(Index = Exogenous_Main),
  freq_metab_main %>% rename(Index = Endogenous_Main)) %>%
  left_join(data_origin %>% select(Index, Origin, Compounds)) %>%
  mutate(Type='Main effect')

freq_age <- dplyr::bind_rows(
  freq_exp_age %>% rename(Index = Exogenous_Age_Interaction),
  freq_metab_age %>% rename(Index = Endogenous_Age_Interaction)) %>%
  left_join(data_origin %>% select(Index, Origin, Compounds)) %>%
  mutate(Type='Metabolite-age interaction')

freq_all <- rbind(freq_main, freq_age) %>%
  select(Index, Compounds, Origin, Type, Frequency, Percentage) %>%
  arrange(desc(Frequency))

write.csv(freq_all, file = '03result/02Table/TabS1.csv', row.names = F)

################ TabS2-3 ################################
load('03result/Ultra_Robust_Biomarkers.RData')
load('01data/02cleandata/data_origin.RData')
data_anno_sel <- data_origin %>%
  filter(Index %in% c(robust_exp_main$Index, robust_metab_main$Index, robust_exp_age$Index, robust_metab_age$Index)) %>%
  select(c("Mode","Index","Compounds",
           "Origin","Class",
           "Formula","Q1..Da.",
           "RT..min.","Adduct",
           "Mass.error","Level","score",
           "PubChem","HMDB_new")) %>%
  rename('Exact m/z' = 'Q1..Da.',
         'RT (min)' = 'RT..min.',
         'Mass Error (ppm)'= 'Mass.error',
         'MSI Level' = 'Level',
         'MS/MS Match Score' = 'score',
         'HMDB' = 'HMDB_new') %>%
  mutate(
    Mode = case_when(
      Mode == "C18_positive" ~ "ESI+",
      Mode == "C18_negative" ~ "ESI-",
      TRUE ~ Mode # 其他情况保持原样
    )
  )

write.csv(data_anno_sel, '03result/02Table/TabS2.csv', row.names = F)

robust_treat_main <- dplyr::bind_rows(
  robust_exp_main,robust_metab_main) %>%
  left_join(data_origin %>% select(Index, Origin, Compounds)) %>%
  mutate(Type = 'Main effect') %>%
  select(Index, Compounds, Origin, Percentage, Type,
         HR_diag, LCI_diag, UCI_diag, P_Value_diag,
         HR_treat,  LCI_treat,  UCI_treat,  P_Value_treat)

robust_treat_age <- dplyr::bind_rows(
  robust_exp_age,robust_metab_age) %>%
  left_join(data_origin %>% select(Index, Origin, Compounds)) %>%
  mutate(Type = 'Metabolite-age interaction') %>%
  select(Index, Compounds, Origin, Percentage, Type,
         HR_diag, LCI_diag, UCI_diag, P_Value_diag,
         HR_treat,  LCI_treat,  UCI_treat,  P_Value_treat)

robust_delay_main <- dplyr::bind_rows(
  robust_exp_main,robust_metab_main) %>%
  left_join(data_origin %>% select(Index, Origin, Compounds)) %>%
  mutate(Type = 'Main effect') %>%
  select(Index, Compounds, Origin, Percentage, Type,
         HR_diag, LCI_diag, UCI_diag, P_Value_diag,
         HR_delay,  LCI_delay,  UCI_delay,  P_Value_delay)

robust_delay_age <- dplyr::bind_rows(
  robust_exp_age, robust_metab_age) %>%
  left_join(data_origin %>% select(Index, Origin, Compounds)) %>%
  mutate(Type = 'Metabolite-age interaction') %>%
  select(Index, Compounds, Origin, Percentage, Type,
         HR_diag, LCI_diag, UCI_diag, P_Value_diag,
         HR_delay,  LCI_delay,  UCI_delay,  P_Value_delay)

robust_treat <- rbind(robust_treat_main, robust_treat_age)
robust_delay <- rbind(robust_delay_main, robust_delay_age)

write.csv(robust_treat, '03result/02Table/TabS3.csv', row.names = F)
write.csv(robust_delay, '03result/02Table/TabS4.csv', row.names = F)




