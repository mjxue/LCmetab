############################### Biomarker scrennning function ###########################
############################### 通用分析函数更新 ########################################
library(survival)
library(splines)

# 1. 主效应 COX 分析函数
maincox_analysis <- function(data_train, data_valid, var_list, covars, time_col = "OS_diag") {
  res <- data.frame()
  # 动态生成生存对象公式
  formula_base <- paste("Surv(", time_col, ", Dead) ~", paste(covars, collapse = " + "))
  
  for (var in var_list) {
    # 训练集分析
    fit_train <- coxph(as.formula(paste(formula_base, "+", var)), data = data_train)
    coef_train <- summary(fit_train)$coefficients
    confint_train <- summary(fit_train)$conf.int
    
    # 验证集分析
    fit_valid <- coxph(as.formula(paste(formula_base, "+", var)), data = data_valid)
    coef_valid <- summary(fit_valid)$coefficients
    confint_valid <- summary(fit_valid)$conf.int
    
    # 提取结果
    res_row <- data.frame(
      Index = var,
      HR_train = coef_train[var, "exp(coef)"],
      CIL_train = confint_train[var, "lower .95"],
      CIU_train = confint_train[var, "upper .95"],
      Z_Value_train = coef_train[var, 'z'],
      P_Value_train = coef_train[var, "Pr(>|z|)"],
      
      HR_valid = coef_valid[var, "exp(coef)"],
      CIL_valid = confint_valid[var, "lower .95"],
      CIU_valid = confint_valid[var, "upper .95"],
      Z_Value_valid = coef_valid[var, 'z'],
      P_Value_valid = coef_valid[var, "Pr(>|z|)"]
    )
    res <- bind_rows(res, res_row)
  }
  
  # 添加校正p值
  res %>%
    mutate(
      fdr_train = p.adjust(P_Value_train, "fdr"),
      bon_train = p.adjust(P_Value_train, "bonferroni")
    )
}


# 2. 交互作用 COX 分析函数
GEinter_analysis <- function(data_train, data_valid, var_list, interaction_var, covars, time_col = "OS_diag") {
  res <- data.frame()
  # 构建基础公式，确保交互项基效应也在模型中
  formula_base <- paste("Surv(", time_col, ", Dead) ~", paste(unique(c(covars, interaction_var)), collapse = " + "))
  
  for (var in var_list) {
    full_formula <- as.formula(paste0(formula_base, " + ", var, " + ", interaction_var, "*", var))
    
    fit_train <- coxph(full_formula, data = data_train)
    coef_train <- summary(fit_train)$coefficients
    confint_train <- summary(fit_train)$conf.int
    
    fit_valid <- coxph(full_formula, data = data_valid)
    coef_valid <- summary(fit_valid)$coefficients
    confint_valid <- summary(fit_valid)$conf.int
    
    interaction_term <- paste0(interaction_var, ":", var)
    row_train <- which(grepl(interaction_term, rownames(coef_train)))
    row_valid <- which(grepl(interaction_term, rownames(coef_valid)))
    
    if(length(row_train) > 0 & length(row_valid) > 0) {
      res_row <- data.frame(
        Index = var,
        Interaction = interaction_var,
        HR_train = coef_train[row_train, "exp(coef)"],
        CIL_train = confint_train[row_train, "lower .95"],
        CIU_train = confint_train[row_train, "upper .95"],
        Z_Value_train = coef_train[row_train, 'z'],
        P_Value_train = coef_train[row_train, "Pr(>|z|)"],
        
        HR_valid = coef_valid[row_valid, "exp(coef)"],
        CIL_valid = confint_valid[row_valid, "lower .95"],
        CIU_valid = confint_valid[row_valid, "upper .95"],
        Z_Value_valid = coef_valid[row_valid, 'z'],
        P_Value_valid = coef_valid[row_valid, "Pr(>|z|)"]
      )
      res <- bind_rows(res, res_row)
    }
  }
  
  res %>%
    mutate(
      fdr_train = p.adjust(P_Value_train, "fdr"),
      bon_train = p.adjust(P_Value_train, "bonferroni")
    )
}

# diag_metab_nonlinear <- nonlinear_screen(metab_train, metab_valid, metab_list, covars_diag, time_col = "OS_diag", event_col = 'Dead', spline_df = 3)
# data_train <- expose_train
# data_valid <- exp_valid
# var_list <- expose_list
# covars <- covars_diag
# time_col = "OS_diag"
# event_col = "Dead"
# spline_df = 3

# 3. 非线性主效应分析函数
nonlinear_screen <- function(data_train, data_valid, var_list, covars, time_col = "OS_diag", event_col = "Dead", spline_df = 3) {
  
  cat(sprintf("开始对 %d 个特征进行非线性样条筛选...\n", length(var_list)))
  results <- list()
  
  # 构建基础协变量字符串
  covar_formula <- paste(covars, collapse = " + ")
  
  for (feature in var_list) {
    # 引入 tryCatch 防止单个极端分布的代谢物导致整个循环崩溃
    tryCatch({
      
      # --- 模型构建 ---
      # 0. 基础模型 (仅含协变量)
      form_base <- as.formula(paste("Surv(", time_col, ", ", event_col, ") ~ ", covar_formula))
      
      # 1. 线性模型 (协变量 + 线性特征)
      form_linear <- as.formula(paste("Surv(", time_col, ", ", event_col, ") ~ ", covar_formula, " + ", feature))
      
      # 2. 非线性模型 (协变量 + 特征的自然样条)
      form_spline <- as.formula(paste("Surv(", time_col, ", ", event_col, ") ~ ", covar_formula, " + ns(", feature, ", df = ", spline_df, ")"))
      
      # --- 训练集评估 (Derivation) ---
      fit_base_train   <- coxph(form_base, data = data_train)
      fit_linear_train <- coxph(form_linear, data = data_train)
      fit_spline_train <- coxph(form_spline, data = data_train)
      
      # 计算 A: 样条特征的整体显著性 P 值 (样条模型 vs 基础模型)
      p_overall_train <- anova(fit_base_train, fit_spline_train)$`Pr(>|Chi|)`[2]
      
      # 计算 B: 纯非线性成分的显著性 P 值 (样条模型 vs 纯线性模型)
      p_nonlinear_train <- anova(fit_linear_train, fit_spline_train)$`Pr(>|Chi|)`[2]
      
      # --- 验证集评估 (Validation) ---
      fit_base_valid   <- coxph(form_base, data = data_valid)
      fit_linear_valid <- coxph(form_linear, data = data_valid)
      fit_spline_valid <- coxph(form_spline, data = data_valid)
      
      # 计算 C: 样条特征在验证集中的整体显著性 P 值
      p_overall_valid <- anova(fit_base_valid, fit_spline_valid)$`Pr(>|Chi|)`[2]
      p_nonlinear_valid <- anova(fit_linear_valid, fit_spline_valid)$`Pr(>|Chi|)`[2]
      
      # 记录结果
      results[[feature]] <- data.frame(
        Index = feature,
        P_Overall_Train = p_overall_train,
        P_NonLinear_Train = p_nonlinear_train,
        P_Overall_Valid = p_overall_valid,
        P_Nonlinear_Valid = p_nonlinear_valid,
        stringsAsFactors = FALSE
      )
      
    }, error = function(e) {
      cat(sprintf("警告: 特征 %s 样条拟合失败，已跳过。\n", feature))
    })
  }
  
  # 合并结果框
  final_res <- do.call(rbind, results)
  rownames(final_res) <- NULL
  
  # 对训练集的整体 P 值进行 FDR (BH) 校正
  final_res$fdr_overall_train <- p.adjust(final_res$P_Overall_Train, method = "BH")
  final_res$fdr_nonLinear_train <- p.adjust(final_res$P_NonLinear_Train, method = "BH")
  
  return(final_res)
}
# ==========================================
# 辅助函数：全队列单因素/协变量调整 COX 分析，提取标准效应量
# ==========================================

# A. 主效应全队列 Cox
full_cox_main <- function(data, feature_list, covars, time_col = "OS_diag", event_col = "Dead") {
  res_list <- lapply(feature_list, function(feat) {
    form <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", feat, " + ", paste(covars, collapse=" + ")))
    fit <- tryCatch(coxph(form, data = data), error = function(e) NULL)
    
    if (is.null(fit)) return(NULL)
    
    sum_fit <- summary(fit)
    coef_mat <- sum_fit$coefficients
    conf_mat <- sum_fit$conf.int
    
    if (feat %in% rownames(coef_mat)) {
      data.frame(
        Index = feat,
        HR = conf_mat[feat, "exp(coef)"],
        LCI = conf_mat[feat, "lower .95"],
        UCI = conf_mat[feat, "upper .95"],
        P_Value = coef_mat[feat, "Pr(>|z|)"]
      )
    } else {
      NULL
    }
  })
  return(do.call(rbind, res_list))
}

# B. 交互作用全队列 Cox
full_cox_inter <- function(data, feature_list, interaction_var, covars, time_col = "OS_diag", event_col = "Dead") {
  res_list <- lapply(feature_list, function(feat) {
    form <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", feat, " * ", interaction_var, " + ", paste(covars, collapse=" + ")))
    fit <- tryCatch(coxph(form, data = data), error = function(e) NULL)
    
    if (is.null(fit)) return(NULL)
    
    sum_fit <- summary(fit)
    coef_mat <- sum_fit$coefficients
    conf_mat <- sum_fit$conf.int
    
    term_names <- rownames(coef_mat)
    inter_term <- term_names[grepl(":", term_names) & grepl(feat, term_names)]
    
    if (length(inter_term) > 0) {
      inter_name <- inter_term[1]
      data.frame(
        Index = feat,
        HR = conf_mat[inter_name, "exp(coef)"],
        LCI = conf_mat[inter_name, "lower .95"],
        UCI = conf_mat[inter_name, "upper .95"],
        P_Value = coef_mat[inter_name, "Pr(>|z|)"]
      )
    } else {
      NULL
    }
  })
  return(do.call(rbind, res_list))
}

# ==========================================
# 辅助函数：批量进行 PH 假设检验 (基于 Schoenfeld 残差)
# ==========================================

# 1. 主效应 PH 检验
check_ph_main <- function(data, feature_list, covars, time_col = "OS_diag", event_col = "Dead") {
  p_values <- sapply(feature_list, function(feat) {
    # 构建 Cox 公式
    form <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", feat, " + ", paste(covars, collapse=" + ")))
    fit <- tryCatch(coxph(form, data = data), error = function(e) NULL)
    
    if(is.null(fit)) return(NA)
    
    # 计算 Schoenfeld 残差检验
    zph <- tryCatch(cox.zph(fit), error = function(e) NULL)
    if(is.null(zph)) return(NA)
    
    # 提取目标特征的 PH P 值
    return(zph$table[feat, "p"])
  })
  
  return(data.frame(Index = feature_list, PH_P_Value = p_values))
}

# 2. 交互效应 PH 检验
check_ph_inter <- function(data, feature_list, interaction_var, covars, time_col = "OS_diag", event_col = "Dead") {
  p_values <- sapply(feature_list, function(feat) {
    # 构建交互 Cox 公式 (注意这里为了严格测试，包含主效应和交互项)
    form <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", feat, " * ", interaction_var, " + ", paste(covars, collapse=" + ")))
    fit <- tryCatch(coxph(form, data = data), error = function(e) NULL)
    
    if(is.null(fit)) return(NA)
    
    zph <- tryCatch(cox.zph(fit), error = function(e) NULL)
    if(is.null(zph)) return(NA)
    
    # 在交互模型中，交互项在 cox.zph 里的名字通常是 "feat:interaction_var" 或者包含冒号
    term_names <- rownames(zph$table)
    inter_term <- term_names[grepl(":", term_names) & grepl(feat, term_names)]
    
    if(length(inter_term) > 0) {
      return(zph$table[inter_term[1], "p"])
    } else {
      return(NA)
    }
  })
  
  return(data.frame(Index = feature_list, PH_P_Value = p_values))
}
