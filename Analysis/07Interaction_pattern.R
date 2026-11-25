####### Interaction pattern##############################
# 1.Plotting the trend chart of interaction effects###########################
rm(list = ls())
library(readxl)
library(dplyr)
library(survival)

setwd('D:/Cancer Model/03LC/01LC_Metab/v3/')
load('01data/02cleandata/metab_compute.RData')
load('01data/02cleandata/data_origin.RData')
load('03result/Ultra_Robust_Biomarkers.RData')

metab_compute$smoke_bin <- ifelse(metab_compute$Smoke==0,0,ifelse(is.na(metab_compute$Smoke),NA,1))

expose_list <- data_origin$Index[data_origin$Origin=='Exogenous']
metab_list <- data_origin$Index[data_origin$Origin=='Endogenous']

covar_list <- c('Age_diagnosis','Age_blood',"Gender","Smoke","Pack_year","Ctype","Cstage",'radio','smoke_bin','bl2di','timediff')
metab_compute <- subset(metab_compute,select = c("ID",covar_list,'OS_diag','OS_bld','Dead','PFS_diag','PFS',data_origin$Index))

####################function###############################
interaction_c2c = function(model, ds, 
                           target_name,
                           moderator_name,
                           xtitle = "Moderator",
                           ytitle = "Effect of Factor",
                           yscale = "b"){
  
  if(is.null(model)) stop("invalid model input, function stopped.")
  library(RColorBrewer)
  library(ggplot2)
  library(ggExtra)
  
  # get positions of variables in ds and model feedback
  target_pos_in_ds = match(target_name, colnames(ds))
  moderator_pos_in_ds = match(moderator_name, colnames(ds))
  
  target_pos_in_model = match(target_name, names(coxrlt$coefficients))
  moderator_pos_in_model = match(moderator_name, names(coxrlt$coefficients))
  interaction_pos_in_model = length(coxrlt$coefficients)
  
  target_independent_coef = coxrlt$coefficients[target_pos_in_model]
  moderator_independent_coef = coxrlt$coefficients[moderator_pos_in_model]
  interaction_coef = coxrlt$coefficients[interaction_pos_in_model]
  
  cov_matrix = coxrlt$var
  
  # adjust b and ci for target according to moderator level
  interacive_modified_bci = function(tb, tv, ib, iv, ticov, mvalue){
    b = tb + ib * mvalue
    se = (tv + iv * mvalue^2 +2 * mvalue * ticov)^0.5
    
    bcl = b - qnorm(0.025, lower.tail=F)*se
    bcu = b + qnorm(0.025, lower.tail=F)*se
    
    rlt = c(mvalue, b, bcl, bcu, se)
    names(rlt) = c("mvalue", "b", "bcl", "bcu", "se")
    
    return(rlt)
  }
  
  
  # generate matrix by a list of moderator values
  bci_matrix = function(vector, 
                        tb=target_independent_coef, 
                        tv=diag(cov_matrix)[target_pos_in_model],
                        ib=interaction_coef,
                        iv=diag(cov_matrix)[interaction_pos_in_model],
                        ticov=cov_matrix[target_pos_in_model, interaction_pos_in_model]){
    bci = lapply(unlist(vector), 
                 function(x){
                   c = interacive_modified_bci(tb, tv, ib, iv, ticov, x)
                   return(c)
                 }
    )
    
    bci= t(data.frame(bci))
    rownames(bci) = NULL
    
    return(as.data.frame(bci))
  }
  
  # generate a list of moderator values for fig
  #moderator_values_for_fig = seq(min(ds[moderator_pos_in_ds]), max(ds[moderator_pos_in_ds]), 0.01)
  moderator_values_for_fig = ds[,which(colnames(ds)==moderator_name)]
  bci_fig = bci_matrix(unlist(moderator_values_for_fig))
  
  bci_fig$significant = (bci_fig[,"bcl"] * bci_fig[,"bcu"]>0)
  
  # adjust effect scale in b or exp(b)
  if(yscale == "exp"){
    bci_fig$effect = exp(bci_fig$b)
    bci_fig$effect_cl = exp(bci_fig$bcl)
    bci_fig$effect_cu = exp(bci_fig$bcu)
  }else{
    bci_fig$effect = bci_fig$b
    bci_fig$effect_cl = bci_fig$bcl
    bci_fig$effect_cu = bci_fig$bcu
  }
  
  find_intercept_base <- function(x_coords, y_coords, y_target = 1) {
    ord <- order(x_coords)
    x_coords <- x_coords[ord]
    y_coords <- y_coords[ord]
    
    cross_indices <- which(diff(sign(y_coords - y_target)) != 0)

    if (length(cross_indices) == 0) {
      return(numeric(0))
    }

    intercepts <- sapply(cross_indices, function(i) {
      x1 <- x_coords[i]
      y1 <- y_coords[i]
      x2 <- x_coords[i+1]
      y2 <- y_coords[i+1]
      x_intercept <- x1 + (x2 - x1) * (y_target - y1) / (y2 - y1)
      return(x_intercept)
    })
    
    return(intercepts)
  }
  
  if(yscale =='exp'){
    x_intercept_cl_base <- find_intercept_base(bci_fig$mvalue, bci_fig$effect_cl, y_target = 1)
    x_intercept_cu_base <- find_intercept_base(bci_fig$mvalue, bci_fig$effect_cu, y_target = 1)
  }else{
    x_intercept_cl_base <- find_intercept_base(bci_fig$mvalue, bci_fig$effect_cl, y_target = 0)
    x_intercept_cu_base <- find_intercept_base(bci_fig$mvalue, bci_fig$effect_cu, y_target = 0)
  }

  all_intercepts_base <- c(x_intercept_cl_base, x_intercept_cu_base)
  
  segment_data <- data.frame(
    x    = all_intercepts_base,
    xend = all_intercepts_base,
    y    = -Inf,               
    yend = max(bci_fig$effect_cu,na.rm = T)        
  )
  
  min_y_value <- floor(min(c(bci_fig$effect_cl, bci_fig$effect_cu), na.rm = TRUE))
  max_y_value <- floor(max(c(bci_fig$effect_cl, bci_fig$effect_cu), na.rm = TRUE))
  
  y_axis_breaks <- seq(from = -1, to = 2, by = 1)
  
  # draw ci band
  ciband <- ggplot()+
    geom_ribbon(aes(x=mvalue, ymin=effect_cl, ymax=effect_cu,fill = ribbon_color),
                data = bci_fig,  fill="#FFD9D8", alpha = 0.6) +
        geom_line(col = "#4874CB",  size = 0.8, data = bci_fig, aes(x=mvalue, y=effect)) + 
    geom_hline(yintercept = ifelse(yscale=="exp", 1, 0), color = "black", size = 0.3, alpha = 1, linetype = "dashed") +
    geom_segment(
      data = segment_data,
      aes(x = x, y = y, xend = xend, yend = yend),
      linetype = "dashed", 
      color = "black",     
      size = 0.3           
    )+
    scale_x_continuous(breaks = c(30,40,50,60,70,80,90), expand = c(0.01,0.01))+
    scale_y_continuous(breaks = y_axis_breaks,expand = c(0.01,0.05), limits = c(-1,2.2))+
    theme_classic() +
    labs(x = 'Age at sampling', y= ytitle) + 
    theme(axis.title = element_text(face = 1,size=10),
          axis.text= element_text(face = "plain",size=8,colour = "black"),
          axis.line = element_line(colour = "black"))
  
  HR_val <- summary(model)$conf.int[interaction_pos_in_model, "exp(coef)"]
  CIL_val <- summary(model)$conf.int[interaction_pos_in_model, "lower .95"]
  CIU_val <- summary(model)$conf.int[interaction_pos_in_model, "upper .95"]
  Pval_raw <- summary(model)$coefficients[interaction_pos_in_model, "Pr(>|z|)"]
  
  beta_str <- sprintf("%.3f", log(HR_val))
  HR_str <- sprintf("%.3f", HR_val)
  
  if(yscale == "exp"){
    CIL_str <- sprintf("%.3f", CIL_val)
    CIU_str <- sprintf("%.3f", CIU_val)
  }else{
    CIL_str <- sprintf("%.3f", log(CIL_val))
    CIU_str <- sprintf("%.3f", log(CIU_val))
  }
  
  Pval_str <- format(Pval_raw, scientific = TRUE, digits = 3)
  
  if(yscale == "exp"){
    label_v <- paste0(
      "paste(italic(HR)[interaction] == ", HR_str, ", ' (95% CI: ", CIL_str, " ~ ", CIU_str, ")')"
    )
    
    label_p <- paste0(
      "paste(italic(P) == '", Pval_str, "')"
    )
    
    
    
  }else{
    label_v <- paste0(
      "paste(italic(beta)[interaction] == ", beta_str, ", ' (95% CI: ", CIL_str, " ~ ", CIU_str, ")')"
    )
    
    label_p <- paste0(
      "paste(italic(P) == '", Pval_str, "')"
    )
  }
  
  ciband1 <- ciband +
    annotate(
      geom = "text",
      x = min(bci_fig$mvalue, na.rm = T),        # <-- ADJUST this x-coordinate
      y = 2.2*1,       # <-- ADJUST this y-coordinate
      label = label_v,
      parse = TRUE,  # <-- This is the crucial part!
      hjust = 0,     # Left-aligns the text
      size = 3     # Adjust size as needed
    )+
    annotate(
      geom = "text",
      x = min(bci_fig$mvalue, na.rm = T),        # <-- ADJUST this x-coordinate
      y = 2.2*0.9,       # <-- ADJUST this y-coordinate
      label = label_p,
      parse = TRUE,  # <-- This is the crucial part!
      hjust = 0,     # Left-aligns the text
      size = 3     # Adjust size as needed
    )
    
  return(ciband1)
  
}


#################### FigS3 #############################
library(patchwork)
library(cowplot)

plots <- list()
ds <- metab_compute
metablist <- c(robust_exp_age$Index, robust_metab_age$Index)
data_origin$Compounds[data_origin$Index=='MW0076752'] <- 'TG(18:2/18:1/20:4)'

for (i in 1:length(metablist)){
  target_name <- metablist[i]
  coxrlt <- coxph(as.formula(paste0('Surv(OS_diag,Dead)~Gender+Smoke+Ctype+Cstage+Age_blood+Age_blood*',target_name)),data = metab_compute)
  moderator_name <- "Age_blood"
  xtitle <- "Age_blood"
  ytitle <- paste0("beta of ",data_origin$Compounds[data_origin$Index==target_name])
  yscale <- "b"
  p = interaction_c2c(coxrlt, ds, 
                      target_name=target_name,
                      moderator_name =moderator_name,
                      xtitle = xtitle,
                      ytitle = ytitle,
                      yscale = yscale)
  plots[[target_name]] <- p
}

final_plot <- plot_grid(
  plotlist = plots, 
  ncol = 3,
  align = "hv",
  labels = "AUTO",            
  label_size = 12,          
  label_fontface = "plain", 
  byrow = FALSE
)

pdf("03result/01Figure/FigS3.pdf",width=9.1, height=3.015)
final_plot
dev.off()

# 2.Cross-table analysis###########################
################ interaction #########################################################
rm(list = ls())
library(readxl)
library(dplyr)
library(survival)

setwd('D:/Cancer Model/03LC/01LC_Metab/v3/')
load('01data/02cleandata/metab_compute.RData')
load('01data/02cleandata/data_origin.RData')
load('03result/Ultra_Robust_Biomarkers.RData')

metab_compute$smoke_bin <- ifelse(metab_compute$Smoke==0,0,ifelse(is.na(metab_compute$Smoke),NA,1))

covar_list <- c('Age_diagnosis','Age_blood',"Gender","Smoke","Pack_year","Ctype","Cstage",'radio','smoke_bin','bl2di','timediff')
metab_compute <- subset(metab_compute,select = c("ID",covar_list,'OS_diag','OS_bld','Dead','PFS_diag','PFS',data_origin$Index))

library(plyr)
library(survival)
ds <- metab_compute

metablist <- c(robust_exp_age$Index, robust_metab_age$Index)

all_results <- list()

for (i in 1:length(metablist)) {
  metab <- metablist[i]
  
  type1_col <- paste0(metab, "_type")
  type2_col <- paste0("Age_type")
  comb_type_col <- paste0("type_", i)
  
  ds[[type1_col]] <- ifelse(ds[[metab]] <= mean(ds[[metab]], na.rm = TRUE), 'Low', 'High')
  ds[[type2_col]] <- ifelse(ds$Age_blood < 65, 'Young', 'Elder')
  
  ds[[type1_col]] <- factor(ds[[type1_col]], levels = c('Low', 'High'))
  ds[[type2_col]] <- factor(ds[[type2_col]], levels = c('Young', 'Elder'))
  
  # 创建组合类型变量
  ds[[comb_type_col]] <- NA
  
  ds[[comb_type_col]][ds[[type1_col]] == 'Low' & ds[[type2_col]] == 'Young'] <- 1
  ds[[comb_type_col]][ds[[type1_col]] == 'Low' & ds[[type2_col]] == 'Elder'] <- 2
  ds[[comb_type_col]][ds[[type1_col]] == 'High' & ds[[type2_col]] == 'Young'] <- 3
  ds[[comb_type_col]][ds[[type1_col]] == 'High' & ds[[type2_col]] == 'Elder'] <- 4


  if (length(which(is.na(ds[[comb_type_col]])=='TRUE'))>0) {
    ds <- ds[-which(is.na(ds[[comb_type_col]])=='TRUE'),]
  }

  ds[[comb_type_col]] <- factor(ds[[comb_type_col]])
  
  formula1 <- as.formula(paste("Surv(OS_diag, Dead) ~ Age_blood + as.numeric(Gender) + as.factor(Smoke) + as.factor(Ctype) + as.factor(Cstage) +", comb_type_col))
  coxrlt1 <- coxph(formula1, data = ds)
  
  HR1 <- c(NA,summary(coxrlt1)$coefficients[, 2])
  CIL1 <- c(NA,summary(coxrlt1)$conf.int[, 3])
  CIU1 <- c(NA,summary(coxrlt1)$conf.int[, 4])
  PValue1 <- c(NA,summary(coxrlt1)$coefficients[, 5])
  
  res1 <- data.frame(HR = HR1, CIL = CIL1, CIU = CIU1, PValue = PValue1)

  type_terms <- grep(paste0("^", comb_type_col), rownames(res1), value = TRUE)
  if (length(type_terms) > 0) {
    res1_selected <- res1[c(1,which(rownames(res1) %in% type_terms)), ]
    # 添加标识信息
    res1_selected$Index <- metab
    res1_selected$Model <- "Combined Type"
    res1_selected$Term <- rownames(res1_selected)

    comb_factor <- factor(ds[[comb_type_col]], levels = 1:4)
    res1_selected$Number <- as.integer(table(comb_factor))
    res1_selected$Dead   <- as.integer(tapply(ds$Dead, comb_factor, sum, default = 0))
    rownames(res1_selected) <- NULL
  } else {
    res1_selected <- data.frame()
  }
  
  formula2 <- as.formula(paste("Surv(OS_diag, Dead) ~ Age_blood + as.numeric(Gender) + as.factor(Smoke) + as.factor(Ctype) + as.factor(Cstage) +", 
                               type1_col, "*", type2_col))
  coxrlt2 <- coxph(formula2, data = ds)
  
  HR2 <- summary(coxrlt2)$coefficients[, 2]
  CIL2 <- summary(coxrlt2)$conf.int[, 3]
  CIU2 <- summary(coxrlt2)$conf.int[, 4]
  PValue2 <- summary(coxrlt2)$coefficients[, 5]
  
  res2 <- data.frame(HR = HR2, CIL = CIL2, CIU = CIU2, PValue = PValue2)
  
  interaction_pattern <- paste0(type1_col, "|", type2_col)
  interaction_terms <- grep(interaction_pattern, rownames(res2), value = TRUE)
  if (length(interaction_terms) > 0) {
    res2_selected <- res2[interaction_terms, ]
    res2_selected$Index <- metab
    res2_selected$Model <- "Interaction"
    res2_selected$Term <- rownames(res2_selected)
    res2_selected$Number <- NA
    res2_selected$Dead <- NA
    rownames(res2_selected) <- NULL
  } else {
    res2_selected <- data.frame()
  }
  
  if (nrow(res1_selected) > 0 & nrow(res2_selected) > 0) {
    combined_res <- rbind(res1_selected, res2_selected)
    combined_res$Pair_Index <- i
    
    combined_res <- combined_res[, c("Index", "Pair_Index", "Model", "Term", 'Number', 'Dead',"HR", "CIL", "CIU", "PValue")]
    
    all_results[[i]] <- combined_res
    
    cat("After the", i, "th analysis was completed, the following was extracted:", nrow(res1_selected), "+", nrow(res2_selected), "=", nrow(combined_res), "\n\n")
  } else {
    cat("The", i, "th analysis encountered an issue and the result was empty\n\n")
  }
  
  ds[[type1_col]] <- NULL
  ds[[type2_col]] <- NULL
  ds[[comb_type_col]] <- NULL
}

final_results <- do.call(rbind, all_results)

final_results_new <- final_results %>%
  left_join(data_origin %>% select('Index','Compounds')) %>%
  mutate( `HR (95% CI)` = ifelse(is.na(HR), 
                                 "Reference", 
                                 paste0(sprintf('%.3f', HR), ' (', sprintf('%.3f', CIL), ', ', sprintf('%.3f', CIU), ')')),
          `Crude mortality` = paste0(sprintf('%.2f',Dead/Number*100),'%'),
          Pair = paste0(Compounds,' × Age')) %>%
  select(Pair, Term, Number, Dead, `Crude mortality`, `HR (95% CI)`, PValue)

write.csv(final_results_new, "03result/02Table/TabS7.csv", row.names = FALSE)

# 3. Age group analysis######################
rm(list = ls())
library(readxl)
library(dplyr)
library(survival)

setwd('D:/Cancer Model/03LC/01LC_Metab/v3/')
source('02code/function/chisq_strata_metab.R')
load('01data/02cleandata/metab_compute.RData')
load('01data/02cleandata/data_origin.RData')
load('03result/Ultra_Robust_Biomarkers.RData')

metab_compute$smoke_bin <- ifelse(metab_compute$Smoke==0,0,ifelse(is.na(metab_compute$Smoke),NA,1))

covar_list <- c('Age_diagnosis','Age_blood',"Gender","Smoke","Pack_year","Ctype","Cstage",'radio','smoke_bin','bl2di','timediff')
metab_compute <- subset(metab_compute,select = c("ID",covar_list,'OS_diag','OS_bld','Dead','PFS_diag','PFS',data_origin$Index))

data_origin$Compounds[data_origin$Index =='MW0076752'] <- 'TG(18:2/18:1/20:4)'

library(plyr)
library(survival)
ds <- metab_compute
ds$UN_group <- ifelse(ds$Age_blood<65,'age<65','age>=65')

metablist <- c(robust_exp_age$Index, robust_metab_age$Index)

all_results <- list()

for (i in 1:length(metablist)) {
  metab <- metablist[i]

  formula1 <- as.formula(paste("Surv(OS_diag, Dead) ~ Age_blood + as.numeric(Gender) + as.factor(Smoke) + as.factor(Ctype) + as.factor(Cstage) +", metab))
  coxrlt1 <- coxph(formula1, data = ds[ds$UN_group=='age<65',])
  coxrlt2 <- coxph(formula1, data = ds[ds$UN_group=='age>=65',])

  res <- dplyr::bind_rows(
    data.frame(Index=metab,Variable='UN standard', Type='age<65',
               beta=summary(coxrlt1)$coefficients[metab, 1],
               se=summary(coxrlt1)$coefficients[metab, 3],
               CIL=log(summary(coxrlt1)$conf.int[metab, 3]),
               CIU=log(summary(coxrlt1)$conf.int[metab, 4]),
               PValue=summary(coxrlt1)$coefficients[metab, 5]),
    
    data.frame(Index=metab,Variable='UN standard', Type='age>=65',
               beta=summary(coxrlt2)$coefficients[metab, 1],
               se=summary(coxrlt2)$coefficients[metab, 3],
               CIL=log(summary(coxrlt2)$conf.int[metab, 3]),
               CIU=log(summary(coxrlt2)$conf.int[metab, 4]),
               PValue=summary(coxrlt2)$coefficients[metab, 5])
  )

  all_results[[i]] <- res
}

final_results <- do.call(rbind, all_results)

final_results_phet <- chisq_strata_metab(data = final_results, metablist = metablist, strata = c('UN standard'))


library(forestplot)
library(dplyr)
library(grid)
library(cowplot)

format_p_plain <- function(p) {
  
  val_str <- sprintf("%.2e", p)
  
  return(parse(text = paste0("plain('", val_str, "')"))[[1]])
}

metablist_unique <- unique(final_results_phet$Index)
forest_list <- list()

for (metab in metablist_unique) {
  
  df_sub <- final_results_phet %>% filter(Index == metab)
  
  labels_col1 <- list("Subgroup")
  labels_col2 <- list("HR (95% CI)")
  labels_col3 <- list("P")
  labels_col4 <- list("Phet")
  
  mean_val <- c(NA)
  lower_val <- c(NA)
  upper_val <- c(NA)
  is_summary <- c(TRUE) 
  
  for (std in c("UN standard")) {
    df_std <- df_sub %>% filter(Variable == std)
    if (nrow(df_std) == 0) next
    
    phet_expr <- sprintf("%.2e", df_std$phet[1])

    labels_col1 <- append(labels_col1, list(std))
    labels_col2 <- append(labels_col2, list(""))
    labels_col3 <- append(labels_col3, list(""))
    labels_col4 <- append(labels_col4, list(phet_expr))
    
    mean_val <- c(mean_val, NA)
    lower_val <- c(lower_val, NA)
    upper_val <- c(upper_val, NA)
    is_summary <- c(is_summary, FALSE) 
    
    for (i in 1:nrow(df_std)) {
      hr <- exp(df_std$beta[i])
      cil <- exp(df_std$CIL[i])
      ciu <- exp(df_std$CIU[i])
      hr_str <- sprintf("%.2f (%.2f - %.2f)", hr, cil, ciu)
      
      pval_expr <- sprintf("%.2e", df_std$PValue[i])
      
      labels_col1 <- append(labels_col1, list(paste0("    ", df_std$Type[i]))) 
      labels_col2 <- append(labels_col2, list(hr_str))
      labels_col3 <- append(labels_col3, list(pval_expr))
      labels_col4 <- append(labels_col4, list("")) 
      
      mean_val <- c(mean_val, hr)
      lower_val <- c(lower_val, cil)
      upper_val <- c(upper_val, ciu)
      is_summary <- c(is_summary, FALSE) 
    }
  }
  
  labeltext <- list(labels_col1, labels_col2, labels_col3, labels_col4)
  compound <- data_origin$Compounds[data_origin$Index == metab]
  
  p <- forestplot(
    labeltext = labeltext,
    mean = mean_val,
    lower = lower_val,
    upper = upper_val,
    is.summary = is_summary,
    zero = 1,                 
    xlog = FALSE,
    clip = c(0.5, 2),
    xticks = c(0.5, 1, 1.5, 2),
    boxsize = 0.25,
    lineheight = unit(0.7, "cm"),
    graphwidth = unit(4, "cm"),  
    colgap = unit(4, "mm"),      
    txt_gp = fpTxtGp(
      label = gpar(cex = 0.8, fontfamily = "sans"),
      ticks = gpar(cex = 0.75),                      
      xlab = gpar(cex = 0.85),
      title = gpar(cex = 0.95, fontface = "bold")
    ),
    
    title = paste(compound), 
    xlab = "Hazard Ratio (95% CI)",
    col = fpColors(box = "#1874CD", line = "#104E8B", zero = "gray50"), 
    lwd.zero = 2,             
    lwd.ci = 1.5,
    graph.pos = 2             
  )
  
  grob_p <- grid::grid.grabExpr(print(p))
  forest_list[[metab]] <- grob_p
}


combined_plot <- plot_grid(
  plotlist = forest_list, 
  nrow = 3,                     
  labels = "AUTO",              
  label_size = 14,
  label_fontface = "plain"
)


ggsave(
  filename = "03result/01Figure/FigS4.pdf", 
  plot = combined_plot, 
  width = 5.83,                   
  height = 6, 
  limitsize = FALSE
)




