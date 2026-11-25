##########################################
# Permutation-based BIF threshold calibration
##########################################

rm(list = ls())
library(readxl)
library(dplyr)
library(survival)
library(foreach)
library(doParallel)  

# ==========================================
# Record start time
# ==========================================
start_time <- Sys.time()
cat(sprintf("========================================================\n"))
cat(sprintf("The program has started to run: %s\n", start_time))
cat(sprintf("========================================================\n"))

setwd('/data1/user/mjxue/02Cancer/01LC/Metabolomics/')
load('01data/data0531/metab_compute.RData')
load('01data/data0531/data_origin.RData')
source('02code/Biomarker_screening_function.R')  

metab_compute$smoke_bin <- ifelse(metab_compute$Smoke==0,0,ifelse(is.na(metab_compute$Smoke),NA,1))

expose_list <- data_origin$Index[data_origin$Origin=='Exogenous']
metab_list <- data_origin$Index[data_origin$Origin=='Endogenous']

covar_list <- c('Age_diagnosis','Age_blood',"Gender","Smoke","Pack_year","Ctype","Cstage",'radio','smoke_bin','bl2di','timediff')
metab_compute <- subset(metab_compute, select = c("ID",covar_list,'OS_diag','OS_bld','Dead','PFS_diag','PFS',data_origin$Index))

covars_diag <- c("Gender", "Smoke", "Ctype", "Cstage", "Age_blood")

# ==========================================
# Parameter settings
# ==========================================
n_cores <- 40
n_bootstrap <- 1000
n_permutations <- 100

# ==========================================
# Step 1: Pre-generate all scrambled indices
# ==========================================

n_total <- nrow(metab_compute)

# Generate all permutations and then scramble the outcome indices
set.seed(1230)
shuffle_indices <- lapply(1:n_permutations, function(p) {
  sample(1:n_total, size = n_total, replace = FALSE)
})

cat(sprintf("========================================================"))
cat(sprintf("Permutation BIF Zero Distribution Construction\n"))
cat(sprintf("Number of randomizations: %d | Each Bootstrap: %d | Parallel cores: %d", n_permutations, n_bootstrap, n_cores))
cat(sprintf("========================================================"))

# ==========================================
# Step 2: Start a cluster on the outer layer
# ==========================================

cl <- makeCluster(n_cores)
registerDoParallel(cl)

clusterExport(cl, c("metab_compute", "expose_list", "metab_list", "covars_diag",
                    "shuffle_indices", "n_total", "n_bootstrap"),
              envir = environment())

clusterEvalQ(cl, {
  library(survival)
  library(dplyr)
})

# ==========================================
# Step 3: Permutation loop
# ==========================================

# Storage zero distribution pool
null_bif_pool <- list(
  exp_main   = numeric(0),
  metab_main = numeric(0),
  exp_age    = numeric(0),
  metab_age  = numeric(0)
)

for (perm_i in 1:n_permutations) {
  loop_start <- Sys.time()
  cat(sprintf("\n--- Permutation %d/%d begeins at %s ---\n", 
              perm_i, n_permutations, format(loop_start, "%Y-%m-%d %H:%M:%S")))
  
  # The currently disrupted index
  shuffle_idx <- shuffle_indices[[perm_i]]
  
  # Parallel Bootstrap: The worker internally reorders the outcomes based on the index, without passing the entire matrix
  results <- foreach(
    i = 1:n_bootstrap, 
    .packages = c("survival", "dplyr"), 
    .export = c("maincox_analysis", "GEinter_analysis")
  ) %dopar% {
    
    set.seed(1230 + perm_i * 10000 + i)     
    
    # ---- Inside the worker: Shuffle the outcome using an index + perform bootstrap sampling ----
    metab_local <- metab_compute
    
    metab_local$OS_diag  <- metab_compute$OS_diag[shuffle_idx]
    metab_local$Dead     <- metab_compute$Dead[shuffle_idx]
    metab_local$OS_bld   <- metab_compute$OS_bld[shuffle_idx]
    metab_local$PFS_diag <- metab_compute$PFS_diag[shuffle_idx]
    metab_local$PFS      <- metab_compute$PFS[shuffle_idx]
    
    boot_index <- sample(1:n_total, size = n_total, replace = TRUE)
    train_data <- metab_local[boot_index, ]
    oob_index <- setdiff(1:n_total, unique(boot_index))
    valid_data <- metab_local[oob_index, ]    
    
    expose_train <- subset(train_data, select = c(expose_list, covars_diag, 'OS_diag', 'Dead'))
    expose_valid <- subset(valid_data, select = c(expose_list, covars_diag, 'OS_diag', 'Dead'))
    metab_train  <- subset(train_data, select = c(metab_list, covars_diag, 'OS_diag', 'Dead'))
    metab_valid  <- subset(valid_data, select = c(metab_list, covars_diag, 'OS_diag', 'Dead'))
    
    diag_exp_main   <- maincox_analysis(expose_train, expose_valid, expose_list, covars_diag, time_col = "OS_diag")
    diag_metab_main <- maincox_analysis(metab_train, metab_valid, metab_list, covars_diag, time_col = "OS_diag")
    
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
    
    diag_exp_age   <- GEinter_analysis(expose_train, expose_valid, expose_list, 
                                       interaction_var = "Age_blood", covars = covars_diag)
    diag_metab_age <- GEinter_analysis(metab_train, metab_valid, metab_list, 
                                       interaction_var = "Age_blood", covars = covars_diag)
    
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
    
    list(
      exp_main_features   = valid_diag_exp_main$Index,
      metab_main_features = valid_diag_metab_main$Index,
      exp_age_features    = valid_diag_exp_age$Index,
      metab_age_features  = valid_diag_metab_age$Index
    )
  }
  
  # ---- Collect this round of BIF and add it to the zero distribution pool ----
  calc_bif_vector <- function(feature_list, all_features) {
    freq <- table(unlist(feature_list))
    bif <- setNames(rep(0, length(all_features)), all_features)
    bif[names(freq)] <- as.numeric(freq) / n_bootstrap * 100
    return(as.numeric(bif))
  }
  
  null_bif_pool$exp_main   <- c(null_bif_pool$exp_main,   calc_bif_vector(lapply(results, function(x) x$exp_main_features), expose_list))
  null_bif_pool$metab_main <- c(null_bif_pool$metab_main, calc_bif_vector(lapply(results, function(x) x$metab_main_features), metab_list))
  null_bif_pool$exp_age    <- c(null_bif_pool$exp_age,    calc_bif_vector(lapply(results, function(x) x$exp_age_features), expose_list))
  null_bif_pool$metab_age  <- c(null_bif_pool$metab_age,  calc_bif_vector(lapply(results, function(x) x$metab_age_features), metab_list))
  
  cat(sprintf("完成. 零分布池大小: exp_main=%d, metab_main=%d, exp_age=%d, metab_age=%d\n",
              length(null_bif_pool$exp_main), 
              length(null_bif_pool$metab_main),
              length(null_bif_pool$exp_age), 
              length(null_bif_pool$metab_age)))
  
  # 计算耗时与剩余时间（你之前已加入，可以保留）
  loop_end <- Sys.time()
  elapsed_this <- difftime(loop_end, loop_start, units = "mins")
  elapsed_total <- difftime(loop_end, start_time, units = "mins")
  avg_per_loop <- elapsed_total / perm_i
  remaining <- avg_per_loop * (n_permutations - perm_i)
  
  cat(sprintf("This round took: %.1f minutes | Average: %.1f 分钟/次", elapsed_this, avg_per_loop))
  cat(sprintf("Remaining time: %.0f minutes (%.1f hours) | Estimated completion: %s", remaining, remaining/60, 
              format(Sys.time() + remaining * 60, "%Y-%m-%d %H:%M:%S")))
}

end_time <- Sys.time()
cat(sprintf("========================================================\n"))
cat(sprintf("All completed! Total time taken: %.1f hours\n", difftime(end_time, start_time, units = "hours")))
cat(sprintf("========================================================\n"))

stopCluster(cl)

# # Save the zero distribution
save(null_bif_pool, file = '03result/permutation_null_boot1000_perm100_1230.RData')


# ==========================================
# Step 4: Determine the unified BIF threshold based on the zero distribution
# ==========================================
load('03result/permutation_null_boot1000_perm100_1230.RData')
cat(sprintf("========================================================"))
cat(sprintf("Zero distribution quantile → Unified BIF threshold"))
cat(sprintf("========================================================"))

quantiles_to_report <- c(0.50, 0.75, 0.90, 0.95, 0.99, 0.999)

calc_thresholds <- function(null_vec, label) {
  cat(sprintf("
--- %s ---
", label))
  cat(sprintf("  The number of zero distribution values: %d
", length(null_vec)))
  cat(sprintf("  Non-zero proportion: %.1f%%
", mean(null_vec > 0) * 100))
  cat(sprintf("  Mean: %.1f%% | Median: %.1f%% | Maximum: %.1f%%
", 
              mean(null_vec), median(null_vec), max(null_vec)))
  
  qvals <- quantile(null_vec, probs = quantiles_to_report)
  for (i in seq_along(quantiles_to_report)) {
    cat(sprintf("  %s quantile: %.1f%%
", 
                paste0(round(quantiles_to_report[i]*100), "%"), qvals[i]))
  }
  return(qvals)
}

thr_exp_main   <- calc_thresholds(null_bif_pool$exp_main,   "Exogenous exposure, main effect")
thr_metab_main <- calc_thresholds(null_bif_pool$metab_main, "Endogenous metabolites, main effect")
thr_exp_age    <- calc_thresholds(null_bif_pool$exp_age,    "Exogenous exposure, metabolite-age interaction")
thr_metab_age  <- calc_thresholds(null_bif_pool$metab_age,  "Endogenous metabolites, metabolite-age interaction")

# ==========================================
# Step 5: Merge zero distribution → Set a unified threshold
# ==========================================

null_all_combined <- c(null_bif_pool$exp_main, 
                       null_bif_pool$metab_main, 
                       null_bif_pool$exp_age, 
                       null_bif_pool$metab_age)

thr_combined <- calc_thresholds(null_all_combined, "Four types combined")

# ==========================================
# # Step 6: Verify whether the four types of zero distributions can be combined
# ==========================================

cat(sprintf("========================================================"))
cat(sprintf("Four types of zero distribution similarity tests (KS test)"))
cat(sprintf("========================================================"))

set.seed(4200)
sample_size <- min(10000, 
                   length(null_bif_pool$exp_main),
                   length(null_bif_pool$metab_main),
                   length(null_bif_pool$exp_age),
                   length(null_bif_pool$metab_age))

s_exp_main   <- sample(null_bif_pool$exp_main, sample_size)
s_metab_main <- sample(null_bif_pool$metab_main, sample_size)
s_exp_age    <- sample(null_bif_pool$exp_age, sample_size)
s_metab_age  <- sample(null_bif_pool$metab_age, sample_size)

ks_tests <- list(
  "exp_main vs metab_main" = ks.test(s_exp_main, s_metab_main)$p.value,
  "exp_main vs exp_age"    = ks.test(s_exp_main, s_exp_age)$p.value,
  "exp_main vs metab_age"  = ks.test(s_exp_main, s_metab_age)$p.value,
  "metab_main vs exp_age"  = ks.test(s_metab_main, s_exp_age)$p.value,
  "metab_main vs metab_age"= ks.test(s_metab_main, s_metab_age)$p.value,
  "exp_age vs metab_age"   = ks.test(s_exp_age, s_metab_age)$p.value
)

for (nm in names(ks_tests)) {
  cat(sprintf("  %s: p = %.4f
", nm, ks_tests[[nm]]))
}

min_ks_p <- min(unlist(ks_tests))
if (min_ks_p > 0.05) {
  cat(sprintf("
The minimum KS p value = %.4f > 0.05, indicating no significant difference among the four types of distributions. It is recommended to use the combined threshold.", min_ks_p))
  use_combined <- TRUE
} else {
  cat(sprintf("
The minimum KS p-value = %.4f ≤ 0.05, indicating that there are differences among the four types of distributions. It is recommended to use separate thresholds for each type.", min_ks_p))
  use_combined <- FALSE
}

# ==========================================
# Step 7: Determine the final threshold (using the maximum value of the zero distribution method)
# ==========================================

cat(sprintf("\n========================================================\n"))
cat(sprintf("The zero distribution maximum value method is used to determine the BIF threshold\n"))
cat(sprintf("========================================================\n"))

get_max_threshold <- function(null_vec, label) {
  max_val <- max(null_vec)
  n_total <- length(null_vec)
  n_nonzero <- sum(null_vec > 0)
  cat(sprintf("  %s:\n", label))
  cat(sprintf("    Total number of zeros: %d | Non-zero numbers: %d (%.1f%%)\n", 
              n_total, n_nonzero, n_nonzero/n_total*100))
  cat(sprintf("    Maximum value: %.1f%%\n", max_val))
  return(max_val)
}

if (use_combined) {
  final_threshold <- get_max_threshold(null_all_combined, "Four types combined")
  final_thresholds <- NULL
  cat(sprintf("\n>>> Global unified BIF threshold: %.1f%% <<<\n", final_threshold))
  cat(sprintf("Meaning: In the %d iterations × across all metabolites, no pure noise BIF exceeded %.1f%%\n",
              n_permutations, final_threshold))
} else {
  final_threshold <- NULL
  final_thresholds <- list(
    exp_main   = get_max_threshold(null_bif_pool$exp_main, "Exogenous main effect"),
    metab_main = get_max_threshold(null_bif_pool$metab_main, "Endogenous main effect"),
    exp_age    = get_max_threshold(null_bif_pool$exp_age, "Exogenous metabolite-age interaction"),
    metab_age  = get_max_threshold(null_bif_pool$metab_age, "Endogenous metabolite-age interaction")
  )
  cat(sprintf("\n>>> Various BIF thresholds <<<\n"))
  cat(sprintf("  Exogenous main effect: %.1f%%\n", final_thresholds$exp_main))
  cat(sprintf("  Endogenous main effect: %.1f%%\n", final_thresholds$metab_main))
  cat(sprintf("  Exogenous metabolite-age interaction:   %.1f%%\n", final_thresholds$exp_age))
  cat(sprintf("  Endogenous metabolite-age interaction:   %.1f%%\n", final_thresholds$metab_age))
}

