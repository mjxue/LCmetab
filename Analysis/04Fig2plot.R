##########################################Analysis######################################
rm(list = ls())
library(readxl)
library(dplyr)
library(survival)
library(ggplot2)
library(forcats)
library(tidyr)
library(ggpattern)

setwd('D:/Cancer Model/03LC/01LC_Metab/v3/')
load('01data/02cleandata/metab_compute.RData')
load('01data/02cleandata/data_origin.RData')
load('03result/freq_summary3100.RData')
load('03result/full_alldata_diag_res_notreat.RData')

source('02code/function/Fig2AB_plot_function.R')

metab_compute$smoke_bin <- ifelse(metab_compute$Smoke==0,0,ifelse(is.na(metab_compute$Smoke),NA,1))

expose_list <- data_origin$Index[data_origin$Origin=='Exogenous']
metab_list <- data_origin$Index[data_origin$Origin=='Endogenous']

covar_list <- c('Age_diagnosis','Age_blood',"Gender","Smoke","Pack_year","Ctype","Cstage",'radio','smoke_bin','bl2di','timediff')
metab_compute <- subset(metab_compute,select = c("ID",covar_list,'OS_diag','OS_bld','Dead','PFS_diag','PFS',data_origin$Index))

################Fig2A##################################################
v_main <- bind_rows(
  freq_exp_main %>% rename(Index=Exogenous_Main),
  freq_metab_main %>% rename(Index=Endogenous_Main)) %>%
  left_join(
    bind_rows(
      full_diag_exp_main,
      full_diag_metab_main)
  )%>%
  mutate(
    Full_Beta = log(HR_diag),
    Effect = "Main effect",
    Label_Name = Index
  )

v_age <- bind_rows(
  freq_exp_age %>% rename(Index=Exogenous_Age_Interaction),
  freq_metab_age %>% rename(Index=Endogenous_Age_Interaction)) %>%
  left_join(
    bind_rows(
      full_diag_exp_age,
      full_diag_metab_age)
  )%>%
  mutate(
    Full_Beta = log(HR_diag),
    Effect = "Metabolite-age interaction",
    Label_Name = Index
  )


df_plot_color <- bind_rows(v_main, v_age) %>%
  mutate(Boot_Percent_Num = as.numeric(gsub("%", "", Percentage))) %>%
  # replace_na(list(Boot_Percent_Num = 0)) %>%
  mutate(Effect = factor(Effect, levels = c("Main effect", "Metabolite-age interaction"))) %>%
  arrange(desc(Boot_Percent_Num)) %>%
  left_join(data_origin %>% select(Index, Compounds, Origin), by = 'Index') %>%
  subset(select=c('Index','Compounds','Origin','Percentage','Boot_Percent_Num','Effect','Full_Beta'))

df_plot_color$Compounds[df_plot_color$Index =='MW0155706'] <- 'Pregnanolone sulphate'
df_plot_color$Compounds[df_plot_color$Index =='MW0138609'] <- 'K-3-(6"Ac)Gal-7-Rha'
df_plot_color$Compounds[df_plot_color$Index =='MW0076752'] <- 'TG(18:2/18:1/20:4)'
df_plot_color$Compounds[df_plot_color$Index =='MW0012126'] <- '11betaOHAn'
df_plot_color$Compounds[df_plot_color$Index =='MW0126302'] <- 'Pyrazole'
df_plot_color$Compounds[df_plot_color$Index =='MW0062277'] <- 'Prostaglandin J2'

df_plot_color$Label_Name <- df_plot_color$Compounds
df_plot_color <- df_plot_color[1:20,]
# =====================================================================
# 3. barplot
# =====================================================================
gold_threshold <- 73.2
max_beta_abs <- max(abs(df_plot_color$Full_Beta), na.rm = TRUE)

pC_gradient <- ggplot(df_plot_color,
                      aes(x = fct_reorder(Label_Name, Boot_Percent_Num,
                                          .desc = TRUE
                      ),
                      y = Boot_Percent_Num,
                      fill = Full_Beta,
                      pattern = Effect            
                      )) +

  geom_col(width = 0.6, alpha = 0.9) +
  geom_col_pattern(
    width = 0.6, alpha = 0.9,
    pattern_fill    = "#8FAADC",      
    pattern_color   = NA,          
    pattern_alpha   = 0.7,      
    pattern_spacing = 0.1,       
    pattern_angle   = 45        
  ) +
  scale_pattern_manual(
    values = c(
      "Main effect"                 = "none",   
      "Metabolite-age interaction"  = "stripe"   
    )
  ) +
  geom_hline(yintercept = gold_threshold, linetype = "dashed", color = "grey30", linewidth = 0.8) +
  annotate("text", x = 18.15, y = gold_threshold + 4, label = "73.2% Threshold",
           color = "grey30", size = 4.5, hjust = 0, fontface = "italic") +

  scale_fill_gradient2(
    low = "#2276B2",
    mid = "white",
    high = "#BB2A33",
    midpoint = 0,
    limits = c(-max_beta_abs, max_beta_abs),
    breaks = c(-0.3, -0.15, 0, 0.15, 0.3)
  ) +

  guides(fill = guide_colorbar(
    title = "Coefficient",
    direction = "horizontal",
    label.position = "top",
    title.vjust = 0.1,
    barheight = unit(0.4, "cm"),   
    barwidth = unit(5, "cm")
  ),
  pattern = guide_legend(
    title = "Effect type",
    direction = "horizontal",
    # title.vjust = 0,
    keyheight = unit(0.4, "cm"),
    keywidth = unit(0.4, "cm"),   
    override.aes = list(pattern_spacing = 0.02)
    )
  ) +

  scale_y_continuous(limits = c(0, 93), expand = c(0, 0), breaks = seq(0, 93, by = 25)) +
  scale_x_discrete(labels = function(x) gsub(" \\(.*\\)", "", x)) +

  labs(
    title = '',
    x = NULL,
    y = "Selection frequency (%)"
  ) +
  geom_text(aes(label = paste0(Boot_Percent_Num),
                y = Boot_Percent_Num + 2),   
            size = 3, vjust = 0, color = "black")+
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position =  c(0.72,0.98),
    legend.box = "horizontal",
    legend.box.just = "bottom",
    legend.background = element_rect(fill = NA),
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(0.5, "cm"),
    legend.title = element_text(size = 9),
    legend.text =  element_text(size = 9),

    axis.text.y = element_text(color = "black", size = 8, margin = margin(r = 5)),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.title =  element_text(color = "black", size = 10),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_blank(),
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.6),
    strip.text = element_text(face = "bold", size = 11, color = "black")
  )

pdf(file = '03result/01Figure/Fig2A.pdf', height = 5, width = 12)
pC_gradient
dev.off()

##############################Fig2B###################################
library(ggplot2)
library(dplyr)
library(ggrepel)
library(ggbeeswarm)
load('03result/Ultra_Robust_Biomarkers.RData')

df_main <- rbind(full_diag_exp_main,full_diag_metab_main)
df_main <- df_main %>%
  mutate(
    logp = -log10(P_Value_diag),
    SigClass = case_when(
      Index %in% c(robust_exp_main$Index, robust_metab_main$Index) ~ "Significant",
      TRUE ~ "Nonsignificant"
    ),
  )

df_age <- rbind(full_diag_exp_age,full_diag_metab_age)
df_age <- df_age %>%
  mutate(
    logp = -log10(P_Value_diag),
    SigClass = case_when(
      Index %in% c(robust_exp_age$Index, robust_metab_age$Index) ~ "Significant",
      TRUE ~ "Nonsignificant"
    ),
  )

df_main <- left_join(df_main, subset(data_origin, select = c('Index','Class','Compounds')), by = "Index") %>% 
  subset(select = c('Index','Compounds','logp','Class','SigClass')) %>%
  dplyr::rename(logp_main = logp, SigClass_main=SigClass)

df_age <- left_join(df_age, subset(data_origin, select = c('Index','Class','Compounds')), by = "Index") %>%
  subset(select = c('Index','Compounds','logp','Class','SigClass')) %>%
  dplyr::rename(logp_age=logp, SigClass_age=SigClass)

df_combined <- left_join(df_main, df_age)

df_combined$Class <- recode(df_combined$Class,
                            "Inorganics & Metallics" = "Inorganics\n& Metallics",
                            "Phenylpropanoids & Polyketides" = "Phenylpropanoids\n& Polyketides"
)


desired_order <- sort(unique(df_combined$Class))

class_info <- data.frame(
  Class = desired_order,
  class_id = 1:length(desired_order)
)


set.seed(123)
metabolite_positions <- df_combined %>%
  distinct(Index, Class) %>%
  left_join(class_info, by = "Class") %>%
  dplyr::group_by(Class) %>%
  dplyr::mutate(x_pos = class_id[1] - 0.5 + runif(n())) %>%
  ungroup() %>%
  select(Index, x_pos)


df_combined_processed <- left_join(df_combined, metabolite_positions, by = "Index")


df_combined_processed$Class <- factor(df_combined_processed$Class, levels = desired_order)

df_combined_final <- df_combined_processed %>%
  mutate(color_group = Class) %>%
  left_join(subset(data_origin,select = c('Index','Origin')))

df_combined_final$Compounds[df_combined_final$Index=='MW0138609'] <- 'K-3-(6"Ac)Gal-7-Rha'
df_combined_final$Compounds[df_combined_final$Index =='MW0062277'] <- 'Prostaglandin J2'
df_combined_final$Compounds[df_combined_final$Index =='MW0012126'] <- '11beta-OHAn'

library(circlize)
library(dplyr)
library(scales)

#####plot#####################
df_plot <- df_combined_final

df_circlize <- df_plot %>%
  group_by(Class) %>%
  dplyr::mutate(x_idx = row_number()) %>% 
  ungroup() %>%
  arrange(Class)

# 强制锁定分类顺序
df_circlize$Class <- factor(df_circlize$Class, levels = unique(df_circlize$Class))
classes <- levels(df_circlize$Class)

colors <- c('#C61E27', '#F79A3A', "#FCC1C8", "#F9C536", 
            "#9FCE80", "#91D1C2FF", "#ADD8E8", "#4DBBD5FF", 
            "#8491B4FF", "#776EAF", "#2C2364")
names(colors) <- classes

sig_threshold <- -log10(0.05)

y_max_main <- max(df_circlize$logp_main, na.rm = TRUE)
y_min_main <- min(df_circlize$logp_main, na.rm = TRUE)
ylim_main <- c(y_max_main * 1.1, y_min_main - 0.05)

y_max_age <- max(df_circlize$logp_age, na.rm = TRUE)
y_min_age <- min(df_circlize$logp_age, na.rm = TRUE)
ylim_age <- c(y_max_age * 1.1, y_min_age - 0.05)


pdf("03result/Fig2B_Manhattan.pdf", width = 6, height = 6)
par(mar = c(0, 0, 0, 0))

circos.clear() 

circos.par(
  "start.degree" = 85, 
  "gap.after" = c(rep(2, length(classes) - 1), 10),
  "cell.padding" = c(0, 0, 0, 0) 
)


circos.initialize(factors = df_circlize$Class, x = df_circlize$x_idx)

circos.track(
  factors = df_circlize$Class, 
  ylim = c(0, 1),           
  track.height = 0.015,      
  bg.border = NA, 
  track.margin = c(0, 0.1), 
  panel.fun = function(x, y) {
    sector.name = CELL_META$sector.index
    xlim = CELL_META$xlim
    
    circos.rect(
      xleft = xlim[1], ybottom = 0,
      xright = xlim[2], ytop = 1,
      col = 'black', 
      border = NA 
    )
    
    circos.text(
      x = CELL_META$xcenter, 
      y = 2, 
      labels = sector.name, 
      facing = "inside", 
      niceFacing = TRUE, 
      adj = c(0.5, 0),      
      cex = 0.5, font = 2
    )
  }
)

r_top <- get.cell.meta.data("cell.top.radius", track.index = 1)
r_bottom <- get.cell.meta.data("cell.bottom.radius", track.index = 1)
draw.sector(
  start.degree = 0, end.degree = 360, 
  rou1 = r_top, rou2 = r_bottom, 
  col = "transparent", border = "black", lwd = 0.8
)


circos.track(
  factors = df_circlize$Class, 
  ylim = ylim_main, 
  bg.border = NA, #'#4D4D4D'
  track.height = 0.22, 
  track.margin = c(0.01, 0.02), 
  panel.fun = function(x, y) {
    sector.name = CELL_META$sector.index
    
    curr_df = df_circlize[df_circlize$Class == sector.name, ]
    valid_df = curr_df[is.finite(curr_df$logp_main), ]
    valid_df = valid_df[order(valid_df$SigClass_main == "Significant"), ]
    
    if(CELL_META$sector.numeric.index == 1) {
      
      y_breaks <- pretty(c(0, y_max_main)) 
      for(br in y_breaks) {
        r_line <- circlize(CELL_META$xcenter, br)[1, "rou"]
        draw.sector(start.degree = 0, end.degree = 360, 
                    rou1 = r_line, rou2 = r_line, 
                    border = "grey90", lwd = 0.5)
      }
    }
    
    if(nrow(valid_df) > 0) {
      pt_shape <- ifelse(valid_df$Origin == "Exogenous", 17, 16)
      
      pt_col <- ifelse(valid_df$SigClass_main == "Significant", "#2D5887", alpha("#9DD9E3", 0.7))
      
      circos.points(
        x = valid_df$x_idx, 
        y = valid_df$logp_main, 
        col = pt_col, 
        pch = pt_shape, 
        cex = 1
      )
      
      sig_df <- valid_df[valid_df$SigClass_main == "Significant", ]
      if(nrow(sig_df) > 0) {
        circos.text(
          x = sig_df$x_idx, 
          y = sig_df$logp_main + (y_max_main * 0.03),
          labels = sig_df$Compounds, 
          facing = "downward",   
          niceFacing = TRUE,     
          adj = c(0, 2),     
          cex = 0.4,             
          col = "#2D5887"     
        )
      }
    }
  }
)

circos.track(
  factors = df_circlize$Class, 
  ylim = ylim_age, 
  bg.border = NA,
  track.height = 0.22, 
  track.margin = c(0.1, 0.05), 
  panel.fun = function(x, y) {
    sector.name = CELL_META$sector.index
    
    curr_df = df_circlize[df_circlize$Class == sector.name, ]
    valid_df = curr_df[is.finite(curr_df$logp_age), ]
    valid_df = valid_df[order(valid_df$SigClass_age == "Significant"), ]
    
    if(CELL_META$sector.numeric.index == 1) {
      
      y_breaks <- pretty(c(0, y_max_age)) 
      for(br in y_breaks) {
        r_line <- circlize(CELL_META$xcenter, br)[1, "rou"]
        draw.sector(start.degree = 0, end.degree = 360, 
                    rou1 = r_line, rou2 = r_line, 
                    border = "grey90", lwd = 0.5)
      }
    }
    
    
    if(nrow(valid_df) > 0) {
      pt_shape <- ifelse(valid_df$Origin == "Exogenous", 17, 16)
      
      pt_col <- ifelse(valid_df$SigClass_age == "Significant", "#2D5887", alpha("#9DD9E3", 0.7))
      
      circos.points(
        x = valid_df$x_idx, 
        y = valid_df$logp_age, 
        col = pt_col, 
        pch = pt_shape, 
        cex = 1
      )
      
      sig_df <- valid_df[valid_df$SigClass_age == "Significant", ]
      if(nrow(sig_df) > 0) {
        circos.text(
          x = sig_df$x_idx, 
          y = sig_df$logp_age + (y_max_age * 0.03), 
          labels = sig_df$Compounds, 
          facing = "downward", 
          niceFacing = TRUE, 
          adj = c(0, 2), 
          cex = 0.4, 
          col = "#2D5887"
        )
      }
    }
  }
)

draw_custom_yaxis <- function(track_idx, max_val, ylim_vals, axis_title) {
  y_breaks <- pretty(c(0, max_val))
  y_breaks <- y_breaks[y_breaks >= 0 & y_breaks <= max_val]
  
  ref_sector <- classes[1]
  ref_x <- get.cell.meta.data("xcenter", sector.index = ref_sector, track.index = track_idx)
  
  breaks_polar <- circlize(rep(ref_x, length(y_breaks)), y_breaks, sector.index = ref_sector, track.index = track_idx)
  bounds_polar <- circlize(c(ref_x, ref_x), ylim_vals, sector.index = ref_sector, track.index = track_idx)
  
  r_breaks <- breaks_polar[, "rou"]
  r_bottom <- bounds_polar[1, "rou"]
  r_top <- bounds_polar[2, "rou"]

  segments(x0 = 0, y0 = r_bottom, x1 = 0, y1 = r_top, col = "black", lwd = 1)
  
  segments(x0 = 0, y0 = r_breaks, x1 = 0.005, y1 = r_breaks, col = "black", lwd = 1)
  
  text(x = 0.022, y = r_breaks, labels = y_breaks, adj = c(0, 0.5), cex = 0.6)
  
  y_mid <- (r_top + r_bottom) / 2  
  text(
    x = -0.03,      
    y = y_mid,       
    labels = axis_title, 
    font = 1, 
    cex = 0.5, 
    srt = 90         
  )
}

draw_custom_yaxis(track_idx = 2, max_val = y_max_main, ylim_vals = ylim_main, axis_title = 'Main effect')
draw_custom_yaxis(track_idx = 3, max_val = y_max_age, ylim_vals = ylim_age, axis_title = 'Metabolite-age\ninteraction')


text(
  x = -0.3, y = 0.03, 
  labels = 'Robust biomarker', 
  font = 2,         
  cex = 0.6,     
  adj = c(0, 0)   
)
legend(
  x = -0.32, y = 0.03, 
  legend = c("Yes", "No"),
  col = c("#2D5887", "#9DD9E3"),
  pch = 16,           
  pt.cex = 1,     
  bty = "n",      
  title = NULL,
  title.font = 1, 
  cex = 0.6,
  # horiz = T
)

text(
  x = 0.05, y = 0.03, 
  labels = "Origin", 
  font = 2,        
  cex = 0.6,     
  adj = c(0, 0) 
)

legend(
  x = 0.03, y = 0.03,  
  legend = c("Exogenous", "Endogenous"),
  col = "grey50",      
  pch = c(17, 16),     
  pt.cex = 1,
  bty = "n",
  title = NULL,      
  cex = 0.6,
  # horiz = T
)

circos.clear()

dev.off()

####################Fig2CD#####################################################
rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

setwd('D:/Cancer Model/03LC/01LC_Metab/v3/')
load('03result/Ultra_Robust_Biomarkers.RData')
load('01data/02cleandata/data_origin.RData')

df_main <- rbind(robust_exp_main, robust_metab_main) %>%
  select(Index,
         # HR_train, LCI_train, UCI_train, P_Value_train,
         # HR_valid, LCI_valid, UCI_valid, P_Value_valid,
         HR_diag, LCI_diag, UCI_diag, P_Value_diag,
         HR_treat,  LCI_treat,  UCI_treat,  P_Value_treat,
         HR_delay,  LCI_delay,  UCI_delay,  P_Value_delay) %>%
  pivot_longer(
    cols = matches("_(diag|treat|delay)$"),   
    names_to = c(".value", "Analysis"),
    names_pattern = "(.*)_(diag|treat|delay)"
  ) %>%
  mutate(Analysis = ifelse(Analysis == "diag", "Primary analysis", 
                           ifelse(Analysis == "treat",'Treatment-adjusted','Delayed-Cox'))) %>%
  left_join(data_origin %>% select(Index, Compounds, Origin), by = 'Index')

df_age <- rbind(robust_exp_age, robust_metab_age) %>%
  select(Index,
         HR_diag, LCI_diag, UCI_diag, P_Value_diag,
         HR_treat,  LCI_treat,  UCI_treat,  P_Value_treat,
         HR_delay,  LCI_delay,  UCI_delay,  P_Value_delay) %>%
  pivot_longer(
    cols = matches("_(diag|treat|delay)$"),   
    names_to = c(".value", "Analysis"),
    names_pattern = "(.*)_(diag|treat|delay)"
  ) %>%
  mutate(Analysis = ifelse(Analysis == "diag", "Primary analysis", 
                           ifelse(Analysis == "treat",'Treatment-adjusted','Delayed-Cox'))) %>%
  left_join(data_origin %>% select(Index, Compounds, Origin), by = 'Index')

df_main$Compounds[df_main$Index == 'MW0138609'] <- 'K-3-(6"Ac)Gal-7-Rha'
df_main$Compounds[df_main$Index =='MW0012126'] <- '11βOHAn'
df_main$Compounds[df_main$Index =='MW0155706'] <- 'Pregnanolone sulphate'
df_age$Compounds[df_age$Index =='MW0076752'] <- 'TG(18:2/18:1/20:4)'

df_main <- df_main %>%
  mutate(
    beta = log(HR),
    CIL = log(LCI),
    CIU = log(UCI),
    direction = ifelse(HR > 1, "Positive", "Negative")
  ) %>%
  arrange(desc(Origin), direction, desc(Compounds))
df_main$Compounds <- factor(df_main$Compounds, levels = unique(df_main$Compounds))

df_age <- df_age %>%
  mutate(
    beta = log(HR),
    CIL = log(LCI),
    CIU = log(UCI),
    direction = ifelse(HR > 1, "Positive", "Negative")
  ) %>%
  arrange(desc(Origin), direction, desc(Compounds))
df_age$Compounds <- factor(df_age$Compounds, levels = unique(df_age$Compounds))

group_colors <- c("Primary analysis" = "#a52a2a", 'Treatment-adjusted' = "#1f78b4" , 
                  'Delayed-Cox' = "#32A838"
)

p_main <- ggplot(df_main, aes(x = HR, y = Compounds, color = Analysis, group = Analysis)) +  # , group = OS
  geom_point(aes(shape = Origin), position = position_dodge(width = 0.7), size = 1.5) +
  geom_errorbar(aes(xmin = LCI, xmax = UCI),
                orientation = "y", 
                position = position_dodge(width = 0.7),
                width = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_color_manual(name = 'Analysis', values = group_colors,
                     breaks = c("Primary analysis", 
                                'Treatment-adjusted',
                                'Delayed-Cox')) +
  # scale_x_continuous(limits = c(-0.7, 0.5), breaks = seq(-0.6, 0.4, 0.2), expand = expansion(mult = c(0, 0))) +
  labs(x = 'Hazard ratio', y = NULL, title = "Main effect") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'bottom', 
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    plot.title = element_text(size = 12, hjust = 0.5)
  )

p_age <- ggplot(df_age, aes(x = HR, y = Compounds, color = Analysis, group = Analysis)) +
  geom_point(aes(shape = Origin), position = position_dodge(width = 0.7), size = 1.5) +
  geom_errorbar(aes(xmin = LCI, xmax = UCI),
                orientation = "y",
                position = position_dodge(width = 0.7), 
                width = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_color_manual(name = 'Analysis', values = group_colors,
                     breaks = c("Primary analysis",
                                'Treatment-adjusted',
                                'Delayed-Cox')) +
  labs(x = 'Hazard ratio', y = NULL, title = "Metabolite-age interaction") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'bottom', 
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    plot.title = element_text(size = 12, hjust = 0.5)
  )


pCD <- p_main + p_age + plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(
    theme = theme(
      legend.position = 'bottom',       
      legend.box = "horizontal",        
      legend.direction = "horizontal",   
      legend.spacing.x = unit(1.5, "cm"),
      plot.margin = margin(0, 0, 0, 0, "cm"),
      legend.margin = margin(0, 0, 0, 0, "cm"),
      legend.spacing = unit(0, "cm"),
      legend.box.margin = margin(0, 0, 0, 0, "cm")
    )
  )

print(pCD)

pdf(file = '03result/01Figure/Fig2CD.pdf', width = 10.08, height = 5)
pCD
dev.off()

