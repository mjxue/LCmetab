##################Fig3A#############################
rm(list = ls())
library(RColorBrewer)
library(corrplot)
library(dplyr)

setwd('D:/Cancer Model/03LC/LC_Metab/v3/')
load('01data/02cleandata/metab_compute.RData')
load('01data/02cleandata/data_origin.RData')
load('03result/Ultra_Robust_Biomarkers.RData')

metab_compute$smoke_bin <- ifelse(metab_compute$Smoke==0, 0, ifelse(is.na(metab_compute$Smoke), NA, 1))

expose_list <- data_origin$Index[data_origin$Origin=='Exogenous']
metab_list  <- data_origin$Index[data_origin$Origin=='Endogenous']

covar_list <- c('Age_diagnosis','Age_blood',"Gender","Smoke","Pack_year","Ctype","Cstage",'radio','bl2di','timediff')
metab_compute <- subset(metab_compute, select = c("ID", covar_list, 'OS_diag','OS_bld','Dead','PFS_diag','PFS', data_origin$Index))

sigmetab_matrix <- dplyr::bind_rows(
  robust_exp_main, robust_exp_age, robust_metab_main, robust_metab_age
) %>%
  left_join(subset(data_origin, select = c(Index, Compounds, Origin, Class))) %>%
  mutate(Compounds = case_when(
    Index == "MW0138609" ~ "K-3-(6\"Ac)\nGal-7-Rha",
    Index == 'MW0012126' ~ '11βOHAn',
    TRUE                 ~ Compounds
  )) %>%
  arrange(Origin, Class, Compounds) %>%
  select(Index, Compounds, Origin, Class)

ds_plots4 <- subset(metab_compute, select = sigmetab_matrix$Index)
colnames(ds_plots4) <- sigmetab_matrix$Compounds

library(ggplot2)
library(tidyr)
library(Hmisc)

cor_res <- rcorr(as.matrix(ds_plots4), type = "pearson")

cor_r <- cor_res$r
cor_p <- cor_res$P

vars <- colnames(ds_plots4)
n <- length(vars)

df_r <- as.data.frame(as.table(cor_r)) %>% 
  dplyr::rename(Value = Freq)

df_p <- as.data.frame(as.table(cor_p)) %>% 
  dplyr::rename(P_val = Freq)

df <- left_join(df_r, df_p, by = c("Var1", "Var2"))

df <- df %>%
  mutate(
    Position = case_when(
      as.numeric(factor(Var1, levels = vars)) < as.numeric(factor(Var2, levels = vars)) ~ "Upper",
      as.numeric(factor(Var1, levels = vars)) > as.numeric(factor(Var2, levels = vars)) ~ "Lower",
      TRUE ~ "Diag"
    ),
    Sig = case_when(
      is.na(P_val) ~ "",
      P_val < 0.001 ~ "***",
      P_val < 0.01 ~ "**",
      P_val < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  filter(Position != "Diag") 

df$Var1 <- factor(df$Var1, levels = rev(vars)) 
df$Var2 <- factor(df$Var2, levels = vars)

p_corr <- ggplot(df, aes(x = Var2, y = Var1)) +
  geom_point(data = filter(df, Position == "Lower"),
             aes(size = abs(Value), fill = Value), 
             shape = 21, color = "white", stroke = 0.5) +
  
  geom_point(data = filter(df, Position == "Upper"),
             aes(size = abs(Value), fill = Value), 
             shape = 21, color = "white", stroke = 0.5) +
  
  geom_text(data = filter(df, Position == "Lower"),
            aes(label = sprintf("%.2f", Value),
                color = ifelse(abs(Value) > 0.6, "white", "black")),
            size = 5) +
  
  geom_text(data = filter(df, Position == "Upper"),
            aes(label = Sig,
                color = ifelse(abs(Value) > 0.6, "white", "black")),
            size = 6, vjust = 0.75) + 
  
  geom_segment(aes(x = 0.5, y = n + 0.5, xend = n + 0.5, yend = 0.5), 
               color = "#BCBCBC", linewidth = 1) +
  
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, 
                       limit = c(-1, 1), name = "Correlation") +
  scale_color_identity() + 
  scale_size_continuous(range = c(6, 16), guide = "none") + 
  
  scale_x_discrete(position = "top") +
  labs(
    caption = "* p < 0.05    ** p < 0.01    *** p < 0.001"
  ) +
  guides(
    fill = guide_colorbar(
      title.hjust = 0,            
      barwidth = unit(5, "cm"),   
      barheight = unit(0.3, "cm") 
    ))+
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "#EBEBEB", linewidth = 1.2), 
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    axis.text.x.top = element_text(angle = 45, hjust = 0, vjust = 0, size = 13, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    
    legend.position = "bottom",
    legend.justification = "left", 

    legend.title = element_text(face = "bold", vjust =1),
    legend.margin = margin(t = 0, b = 0),
    
    plot.caption = element_text(size = 11, face = "italic", color = "black", 
                                hjust = 1, margin = margin(t = -23)),

    plot.margin = margin(t = 10, r = 10, b = 20, l = 10),
    # ----------------------------------------------
    
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

print(p_corr)

ggsave("03result/01Figure/Fig3A.pdf", plot = p_corr, width = 8, height = 8.17, units = "in", dpi = 300)

