
rm(list = ls())
setwd('D:/Cancer Model/03LC/LC_Metab/v3/')

library(dplyr)
library(tidyr)
library(survival)
library(survminer)
library(CMAverse)

load('03result/Ultra_Robust_Biomarkers.RData')
load('01data/02cleandata/data_origin.RData')
load('01data/02cleandata/metab_compute.RData')

# 基础数据清洗
metab_compute$smoke_bin <- ifelse(metab_compute$Smoke==0, 0, ifelse(is.na(metab_compute$Smoke), NA, 1))
covar_list <- c('Age_diagnosis','Age_blood',"Gender","Smoke","Pack_year","Ctype","Cstage",'radio','smoke_bin','bl2di','timediff')
metab_compute <- subset(metab_compute, select = c("ID", covar_list, 'OS_diag','OS_bld','Dead','PFS_diag','PFS', data_origin$Index))
metab_compute$Ctype <- ifelse(metab_compute$Ctype=='LUAD',0,
                              ifelse(metab_compute$Ctype=='LUSC',1,
                                     ifelse(metab_compute$Ctype=='Other',2,
                                            ifelse(metab_compute$Ctype=='Other_NSCLC',3,NA))))
metab_compute$Ctype <- as.factor(metab_compute$Ctype)

exposelist <- c(robust_exp_main$Index, robust_exp_age$Index)
metablist <- c(robust_metab_main$Index, robust_metab_age$Index)

results_list <- list()
covars_str <- "Age_blood + Gender + Smoke + Ctype + Cstage"
age_cut <- 65

set.seed(2026)

for (i in 1:length(exposelist)) {
  expose <- exposelist[i]
  
  for (j in 1:length(metablist)) {
    metab <- metablist[j]
    
    df <- subset(metab_compute, select = c('OS_diag', 'Dead', 'Age_blood', 'Gender', 'Smoke', 'Ctype', 'Cstage', expose, metab))
    df <- na.omit(df)
    
    need_strata <- (expose %in% robust_exp_age$Index) | (metab %in% robust_metab_age$Index)
    
    if (!need_strata) {
      mx_form <- paste0(metab, " ~ ", expose, " + ", covars_str)
      fit_mx <- lm(as.formula(mx_form), data = df)
      
      ym_form <- paste0("Surv(OS_diag, Dead) ~ ", metab, " + ", covars_str)
      fit_ym <- coxph(as.formula(ym_form), data = df)
      
      res <- tryCatch({
        cmest(
          data = df, model = "gformula", full = TRUE,
          outcome = "OS_diag", event = "Dead",
          exposure = expose, mediator = metab,
          basec = c("Age_blood", "Gender", "Smoke", "Ctype", "Cstage"),
          EMint = FALSE,
          mreg = list("linear"), yreg = "coxph",
          astar = 0, a = 1, mval = list(0),
          estimation = "imputation", inference = "bootstrap",
          nboot = 1000, boot.ci.type = 'per'
        )
      }, error = function(e) {
        message(paste("Overall error:", expose, "-", metab, ":", e$message))
        return(NULL)
      })
      
      if (!is.null(res)) {
        sum_res <- summary(res)
        results_list[[length(results_list) + 1]] <- data.frame(
          
          exposure = expose, mediator = metab, Subgroup = "Overall",
          mx_Est = summary(fit_mx)$coefficients[expose,'Estimate'],
          mx_CIL = confint(fit_mx)[expose,'2.5 %'],
          mx_CIU = confint(fit_mx)[expose,'97.5 %'],
          mx_P = summary(fit_mx)$coefficients[expose,'Pr(>|t|)'],
          ym_Est = summary(fit_ym)$coefficients[metab,'coef'],
          ym_CIL = log(summary(fit_ym)$conf.int[metab,'lower .95']),
          ym_CIU = log(summary(fit_ym)$conf.int[metab,'upper .95']),
          ym_P = summary(fit_ym)$coefficients[metab,'Pr(>|z|)'],
          Rte_Est = sum_res$summarydf['Rte','Estimate'],
          Rte_CIL = sum_res$summarydf['Rte','95% CIL'],
          Rte_CIU = sum_res$summarydf['Rte','95% CIU'],
          Rte_P   = sum_res$summarydf['Rte','P.val'],
          Rpnde_Est = sum_res$summarydf['Rpnde','Estimate'],
          Rpnde_CIL = sum_res$summarydf['Rpnde','95% CIL'],
          Rpnde_CIU = sum_res$summarydf['Rpnde','95% CIU'],
          Rpnde_P   = sum_res$summarydf['Rpnde','P.val'],
          Rpnie_Est = sum_res$summarydf['Rpnie','Estimate'],
          Rpnie_CIL = sum_res$summarydf['Rpnie','95% CIL'],
          Rpnie_CIU = sum_res$summarydf['Rpnie','95% CIU'],
          Rpnie_P   = sum_res$summarydf['Rpnie','P.val'],
          pm_Est = sum_res$summarydf['pm','Estimate'],
          pm_CIL = sum_res$summarydf['pm','95% CIL'],
          pm_CIU = sum_res$summarydf['pm','95% CIU'],
          pm_P   = sum_res$summarydf['pm','P.val'],
          stringsAsFactors = FALSE
        )
      }
      
    } else {

      df$age_group <- ifelse(df$Age_blood >= age_cut, "Age≥65", "Age<65")
      
      for (group in c("Age<65", "Age≥65")) {
        sub_df <- subset(df, age_group == group)
        if (nrow(sub_df) < 50) next   
        
        mx_form <- paste0(metab, " ~ ", expose, " + ", covars_str)
        fit_mx <- lm(as.formula(mx_form), data = sub_df)
        
        ym_form <- paste0("Surv(OS_diag, Dead) ~ ", metab, " + ", covars_str)
        fit_ym <- coxph(as.formula(ym_form), data = sub_df)
        
        res_sub <- tryCatch({
          cmest(
            data = sub_df, model = "gformula", full = TRUE,
            outcome = "OS_diag", event = "Dead",
            exposure = expose, mediator = metab,
            basec = c("Age_blood", "Gender", "Smoke", "Ctype", "Cstage"),
            EMint = FALSE,
            mreg = list("linear"), yreg = "coxph",
            astar = 0, a = 1, mval = list(0),
            estimation = "imputation", inference = "bootstrap",
            nboot = 1000, boot.ci.type = 'per'
          )
        }, error = function(e) {
          message(paste("Strata error:", expose, "-", metab, group, ":", e$message))
          return(NULL)
        })
        
        if (!is.null(res_sub)) {
          sum_sub <- summary(res_sub)
          results_list[[length(results_list) + 1]] <- data.frame(
            exposure = expose, mediator = metab, Subgroup = group,
            mx_Est = summary(fit_mx)$coefficients[expose,'Estimate'],
            mx_CIL = confint(fit_mx)[expose,'2.5 %'],
            mx_CIU = confint(fit_mx)[expose,'97.5 %'],
            mx_P = summary(fit_mx)$coefficients[expose,'Pr(>|t|)'],
            ym_Est = summary(fit_ym)$coefficients[metab,'coef'],
            ym_CIL = log(summary(fit_ym)$conf.int[metab,'lower .95']),
            ym_CIU = log(summary(fit_ym)$conf.int[metab,'upper .95']),
            ym_P = summary(fit_ym)$coefficients[metab,'Pr(>|z|)'],
            Rte_Est = sum_sub$summarydf['Rte','Estimate'],
            Rte_CIL = sum_sub$summarydf['Rte','95% CIL'],
            Rte_CIU = sum_sub$summarydf['Rte','95% CIU'],
            Rte_P   = sum_sub$summarydf['Rte','P.val'],
            Rpnde_Est = sum_sub$summarydf['Rpnde','Estimate'],
            Rpnde_CIL = sum_sub$summarydf['Rpnde','95% CIL'],
            Rpnde_CIU = sum_sub$summarydf['Rpnde','95% CIU'],
            Rpnde_P   = sum_sub$summarydf['Rpnde','P.val'],
            Rpnie_Est = sum_sub$summarydf['Rpnie','Estimate'],
            Rpnie_CIL = sum_sub$summarydf['Rpnie','95% CIL'],
            Rpnie_CIU = sum_sub$summarydf['Rpnie','95% CIU'],
            Rpnie_P   = sum_sub$summarydf['Rpnie','P.val'],
            pm_Est = sum_sub$summarydf['pm','Estimate'],
            pm_CIL = sum_sub$summarydf['pm','95% CIL'],
            pm_CIU = sum_sub$summarydf['pm','95% CIU'],
            pm_P   = sum_sub$summarydf['pm','P.val'],
            stringsAsFactors = FALSE
          )
        }
      }
    }
    
    print(paste("Completed:", expose, "->", metab))
  }
}


final_results <- do.call(rbind, results_list)


media_results <- final_results %>%
  dplyr::rename(Index = exposure) %>% 
  left_join(subset(data_origin, select=c(Index,Compounds)), by = "Index") %>%
  
  dplyr::rename(Exposure = Index, Exposure_name = Compounds) %>%
  dplyr::rename(Index = mediator) %>% 
  left_join(subset(data_origin, select=c(Index,Compounds)), by = "Index") %>%
 
  dplyr::rename(Mediator = Index, Mediator_name = Compounds) %>%
  mutate(Exposure_name = ifelse(Exposure == 'MW0138609', 
                                'K-3-(6"Ac)Gal-7-Rha', 
                                Exposure_name),
         Mediator_name = case_when(
           Mediator == 'MW0155706' ~ 'Pregnanolone sulphate',
           Mediator == 'MW0076752' ~ 'TG(18:2/18:1/20:4)',
           Mediator == 'MW0012126' ~ '11betaOHAn',
           Mediator == 'MW0062277' ~ 'Prostaglandin J2',
           TRUE ~ Mediator_name),
         Rpnie_FDR = p.adjust(Rpnie_P, method = "fdr"),
         Rte_Est=log(Rte_Est), 
         Rte_CIL=log(Rte_CIL), 
         Rte_CIU=log(Rte_CIU), 
         Rpnde_Est=log(Rpnde_Est), 
         Rpnde_CIL=log(Rpnde_CIL), 
         Rpnde_CIU=log(Rpnde_CIU),
         Rpnie_Est=log(Rpnie_Est), 
         Rpnie_CIL=log(Rpnie_CIL), 
         Rpnie_CIU=log(Rpnie_CIU)
  ) %>%
  select(Exposure, Exposure_name, Mediator, Mediator_name, Subgroup,
         mx_Est, mx_CIL, mx_CIU, mx_P, 
         ym_Est, ym_CIL, ym_CIU, ym_P, 
         Rte_Est, Rte_CIL, Rte_CIU, Rte_P, 
         Rpnde_Est, Rpnde_CIL, Rpnde_CIU, Rpnde_P,
         Rpnie_Est, Rpnie_CIL, Rpnie_CIU, Rpnie_P,
         pm_Est, pm_CIL, pm_CIU, pm_P,
         Rpnie_FDR)

media_results$Subgroup[media_results$Subgroup=='Age≥65'] <- 'Age>=65'
media_sig <- media_results[media_results$Rpnie_FDR < 0.05,]

write.csv(media_results, '03result/02Table/TabS8.csv', row.names = F)


####################### Forestplot ###################################################
rm(list = ls())
load('03result/media_results_strata_nointeract.RData')
load('03result/Ultra_Robust_Biomarkers.RData')
load('01data/02cleandata/data_origin.RData')

library(circlize)
library(RColorBrewer)
library(dplyr)
library(forestploter)


df_plot <- media_results %>% filter(Subgroup %in% c('Age>=65','Age<65')) %>%
  mutate(
    Category = Mediator_name,        
    Group = Exposure_name,           
    Value = pm_Est * 100,       
    significance = case_when(
      Rpnie_FDR < 0.001 ~ "***",
      Rpnie_FDR < 0.01 ~ "**",
      Rpnie_FDR < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  arrange(Category, Group) %>%
  group_by(Exposure_name, Mediator_name) %>%
  filter(any(Rpnie_FDR <= 0.05)) %>%
  ungroup()

df_plot_new <- df_plot %>%
  select(Exposure_name, Mediator_name, Subgroup, Rpnie_Est, Rpnie_CIL, Rpnie_CIU, Rpnie_P) %>%
  rename(Exposure = Exposure_name, Mediator = Mediator_name) %>%
  arrange(Exposure, Mediator) %>%
  mutate(
    forest = paste(rep(" ", 20), collapse = " "),
    NIE = as.numeric(sprintf('%.3f',Rpnie_Est)),
    CIL = as.numeric(sprintf('%.3f',Rpnie_CIL)),
    CIU = as.numeric(sprintf('%.3f',Rpnie_CIU)),
    `NIE (95%CI)` = paste0(NIE,' (',CIL,' ~ ',CIU,')'),
    P = Rpnie_P,
    across(
      .cols = c(Exposure, Mediator), 
      .fns = ~ifelse(row_number() %% 2 == 0, '', .)
      )
    ) %>%
  select(Exposure, Mediator, Subgroup, forest, NIE, CIL, CIU, `NIE (95%CI)`, P)



tm <- forest_theme(base_size = 4,
                   refline_lty = "solid",
                   ci_pch = c(15),
                   ci_col = c("#262626"),
                   ci_fill = c('#92D050'),
                   ci_Theight = 0.2,
                   footnote_col = "blue",
                   core = list(
                     padding = unit(c(6, 1), "mm"),
                     bg_params = list(
                       fill = c("#ffffff", "#ffffff") 
                       )
                     ),
                   heading = list(
                     padding = unit(c(6, 1), "mm")  
                     )
                   )


n_rows <- nrow(df_plot_new)
pD_age <- forest(df_plot_new[,c(1, 2, 3, 4, 8, 9)],
                  est = list(df_plot_new$NIE
                  ),
                  lower = list(df_plot_new$CIL
                  ), 
                  upper = list(df_plot_new$CIU
                  ),
                  ci_column = c(4), 
                  ref_line = 0,
                  xlim = c(-0.3,0.2),
                  ticks_at = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2), 
                  nudge_y = 0,
                  theme = tm,
                 row_heights = rep(unit(2, "mm"), n_rows),
                 col_widths =  c(120, 120, 80, 160, 120, 60))


pD_age <- add_border(pD_age, part = "header", where = 'top', gp = gpar(lwd = 1, col = "black"))
pD_age <- add_border(pD_age, part = "header", where = 'bottom', gp = gpar(lwd = 0.75, col = "black"))
pD_age <- add_border(pD_age, part = "body",row = 16, where = 'bottom', gp = gpar(lwd = 1, col = "black"))


pdf("03result/01Figure/media_age_sig_plot_Fig3E.pdf", width = 4.1, height = 1.8)
print(pD_age)
dev.off()



