rm(list = ls())

library(dplyr)
library(survival)
library(rms)
library(ggplot2)
library(patchwork)

# ============================================================
# 0. load data
# ============================================================
setwd('D:/Cancer Model/03LC/01LC_Metab/v3/')
load('01data/02cleandata/metab_compute.RData')
load('01data/02cleandata/data_origin.RData')
load('03result/Ultra_Robust_Biomarkers.RData')

metablist_main <- c(robust_exp_main$Index, robust_metab_main$Index)

covar_list <- c('Age_blood', "Gender", "Smoke", "Ctype", "Cstage")

ds <- metab_compute %>%
  dplyr::select(ID, all_of(covar_list), OS_diag, Dead, all_of(metablist_main)) %>%
  na.omit()   

ds$Gender <- factor(ds$Gender)
ds$Smoke  <- factor(ds$Smoke)
ds$Ctype  <- factor(ds$Ctype)
ds$Cstage <- factor(ds$Cstage)

dd <- datadist(ds)
options(datadist = "dd")

data_origin$Compounds[data_origin$Index == 'MW0138609'] <- 'K-3-(6"Ac)Gal-7-Rha'
data_origin$Compounds[data_origin$Index =='MW0012126'] <- '11βOHAn'
# ============================================================
# 1. RCS analysis function
# ============================================================

rcs_single_metab <- function(data, metab, covars, nk = 3) {

  covar_str <- paste(covars, collapse = " + ")
  formula_str <- paste0("Surv(OS_diag, Dead) ~ rcs(", metab, ", ", nk, ") + ", covar_str)
  f <- as.formula(formula_str)
  
  model <- cph(f, data = data, x = TRUE, y = TRUE, surv = TRUE)
  
  anova_res <- anova(model)
  
  p_col <- which(colnames(anova_res) %in% c("P", "Pr(>|Chi|)"))
  if (length(p_col) == 0) {
    p_col <- ncol(anova_res)
  } else {
    p_col <- p_col[1]
  }
  
  P_overall <- as.numeric(anova_res[1, p_col])
  if(length(P_overall) == 0) P_overall <- NA 
  
  if ("Nonlinear" %in% rownames(anova_res)) {
    P_nonlinear <- as.numeric(anova_res["Nonlinear", p_col])
  } else if (nrow(anova_res) >= 2) {
    P_nonlinear <- as.numeric(anova_res[2, p_col])
  } else {
    P_nonlinear <- NA
  }
  if(length(P_nonlinear) == 0) P_nonlinear <- NA
  
  pred <- Predict(model, name = metab, ref.zero = TRUE, fun = identity)
  pred_df <- as.data.frame(pred)
  names(pred_df)[1] <- "x" 

  compound_name <- data_origin$Compounds[data_origin$Index == metab]
  if (length(compound_name) == 0 || is.na(compound_name[1])) {
    compound_name <- metab
  } else {
    compound_name <- as.character(compound_name[1])
  }
  
  list(
    metab = as.character(metab),
    compound = compound_name,
    model = model,
    P_overall = P_overall,
    P_nonlinear = P_nonlinear,
    pred_df = pred_df,
    n_events = as.integer(model$n['Event']),
    n_total = as.integer(sum(model$n)),
    ref_val = median(data[[metab]], na.rm = TRUE)
  )
}

# ============================================================
# 2. run RCS analysis
# ============================================================
rcs_results <- list()

for (i in seq_along(metablist_main)) {
  metab <- metablist_main[i]
  cat(sprintf("\n========== analysis %d/%d: %s ==========\n", i, length(metablist_main), metab))
  
  res <- rcs_single_metab(
    data = ds,
    metab = metab,
    covars = covar_list,
    nk = 3
  )
  
  cat(sprintf("  Compound: %s\n", res$compound))
  cat(sprintf("  P_overall  : %.4e\n", res$P_overall))
  cat(sprintf("  P_nonlinear: %.4e\n", res$P_nonlinear))
  cat(sprintf("  Events: %d / %d\n", res$n_events, res$n_total))
  
  rcs_results[[i]] <- res
}

# ============================================================
# 3. results
# ============================================================
summary_table <- data.frame(
  Compound    = vapply(rcs_results, function(x) x$compound, character(1)),
  Index       = vapply(rcs_results, function(x) x$metab, character(1)),
  P_overall   = vapply(rcs_results, function(x) x$P_overall, numeric(1)),
  P_nonlinear = vapply(rcs_results, function(x) x$P_nonlinear, numeric(1)),
  N_events    = vapply(rcs_results, function(x) x$n_events, numeric(1)),
  N_total     = vapply(rcs_results, function(x) x$n_total, numeric(1)),
  stringsAsFactors = FALSE
)

print(summary_table)

summary_table$Nonlinear_significant <- ifelse(summary_table$P_nonlinear < 0.05, "Yes", "No")

# write.csv(summary_table, "03result/02Table/TabSx_RCS_summary.csv", row.names = FALSE)

# ============================================================
# 4. plot RCS
# ============================================================

plot_rcs_overlay <- function(rcs_res, raw_data) {
  
  df <- rcs_res$pred_df
  metab <- rcs_res$metab
  p_overall   <- rcs_res$P_overall
  p_nonlinear <- rcs_res$P_nonlinear
  compound    <- rcs_res$compound
  
  df_raw <- data.frame(x = raw_data[[metab]])
  x_limits <- range(df$x, na.rm = TRUE)
  df_raw <- subset(df_raw, x >= x_limits[1] & x <= x_limits[2])
  
  y_min <- min(df$lower, df$yhat, na.rm = TRUE)
  y_max <- max(df$upper, df$yhat, na.rm = TRUE)
  y_span <- y_max - y_min
  
  hist_ratio <- 0.8 
  hist_height <- y_span * hist_ratio
  
  hist_info <- hist(df_raw$x, breaks = 30, plot = FALSE)
  hist_df <- data.frame(
    xmin = hist_info$breaks[-length(hist_info$breaks)],
    xmax = hist_info$breaks[-1],
    count = hist_info$counts
  )
  max_count <- max(hist_df$count, na.rm = TRUE)
  if (max_count == 0) max_count <- 1
  
  hist_df$ymin <- y_min
  hist_df$ymax <- y_min + (hist_df$count / max_count) * hist_height

  format_p <- function(p) {
    if (is.na(p) || is.nan(p)) return("NA")
    sprintf("%.2E", p)
  }
  label_text <- sprintf("P for overall = %s\nP for nonlinear = %s",
                        format_p(p_overall), format_p(p_nonlinear))
  
  x_pos <- x_limits[1] + 0.05 * (x_limits[2] - x_limits[1])
  y_pos <- y_max - 0.05 * y_span
  
  p <- ggplot() +
    geom_rect(
      data = hist_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "grey75", color = "white", linewidth = 0.2, alpha = 0.6
    ) +
    
    geom_ribbon(
      data = df,
      aes(x = x, ymin = lower, ymax = upper),
      fill = "#FFD9D8", alpha = 0.5
    ) +
    
    geom_line(
      data = df,
      aes(x = x, y = yhat),
      color = "#4874CB", linewidth = 1.0
    ) +
    
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    
    annotate("text", x = x_pos, y = y_pos, 
             label = label_text, hjust = 0, vjust = 1, size = 3.5, color = "black") +
    
    labs(title = compound, x = "Metabolite concentration") +
    scale_x_continuous(expand = c(0, 0)) +
    
    scale_y_continuous(
      name = "ln(Hazard Ratio)",
      expand = expansion(mult = c(0, 0.05)),
      sec.axis = sec_axis(
        trans = ~ (. - y_min) / hist_height * max_count,
        name = "Frequency",
        breaks = scales::pretty_breaks(n = 4)
      )
    ) +
    coord_cartesian(xlim = x_limits, ylim = c(y_min, y_max)) +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
      axis.title.x = element_text(size = 10),
      axis.title.y.left = element_text(size = 10, color = "black"),
      axis.title.y.right = element_text(size = 10, color = "black"),
      axis.text.x  = element_text(size = 8, color = "black"),
      axis.text.y.left = element_text(size = 8, color = "black"),
      axis.text.y.right = element_text(size = 8, color = "black"),
      plot.title   = element_text(size = 11, face = "bold", hjust = 0.5),
      plot.margin  = margin(10, 10, 10, 10)
    )
  
  return(p)
}

plot_list <- lapply(rcs_results, function(res) plot_rcs_overlay(res, raw_data = ds))

combined_plot <- wrap_plots(plot_list, ncol = 3, nrow = 2) +
  plot_annotation(
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

print(combined_plot)

ggsave("03result/01Figure/FigS2.pdf", combined_plot, width = 11, height = 7, dpi = 300)
