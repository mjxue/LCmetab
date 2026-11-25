rm(list = ls())
library(readxl)
library(dplyr)
library(survival)
library(tidyr)
library(bnlearn)
library(Rgraphviz)

setwd('D:/Cancer Model/03LC/LC_Metab/v3/')
load('01data/02cleandata/metab_compute.RData')
load('01data/02cleandata/data_origin.RData')
load('03result/Ultra_Robust_Biomarkers.RData')

metab_compute$smoke_bin <- ifelse(metab_compute$Smoke == 0, 0,
                                  ifelse(is.na(metab_compute$Smoke), NA, 1))

covar_list <- c('Age_diagnosis', 'Age_blood', "Gender", "Smoke", "Pack_year",
                "Ctype", "Cstage", 'radio', 'smoke_bin', 'bl2di', 'timediff')
metab_compute <- subset(metab_compute,
                        select = c("ID", covar_list, 'OS_diag', 'OS_bld', 'Dead', 'PFS_diag', 'PFS',
                                   data_origin$Index))

explist   <- c(robust_exp_main$Index, robust_exp_age$Index)
metablist <- c(robust_metab_main$Index, robust_metab_age$Index)
selected_metabs <- c(explist, metablist)

clinical_vars <- c("Age_blood", "Gender", "Smoke", "Cstage")
adjust_vars   <- c("Age_blood", "Gender", "Smoke", "Cstage")
surv_nodes    <- c("Surv_1yr", "Surv_3yr", "Surv_5yr")

net_data <- metab_compute %>%
  select(all_of(clinical_vars), all_of(selected_metabs), OS_diag, Dead)

net_data <- net_data %>% drop_na(all_of(clinical_vars), all_of(selected_metabs))


net_data <- net_data %>%
  mutate(
    os_year = OS_diag,
    Surv_1yr = case_when(
      os_year >= 1                    ~ "Yes",
      os_year < 1 & Dead == 0         ~ "Yes",
      os_year < 1 & Dead == 1         ~ "No",
      TRUE                            ~ NA_character_
    ),
    Surv_3yr = case_when(
      os_year >= 3                    ~ "Yes",
      os_year < 3 & Dead == 0         ~ "Yes",
      os_year < 3 & Dead == 1         ~ "No",
      TRUE                            ~ NA_character_
    ),
    Surv_5yr = case_when(
      os_year >= 5                    ~ "Yes",
      os_year < 5 & Dead == 0         ~ "Yes",
      os_year < 5 & Dead == 1         ~ "No",
      TRUE                            ~ NA_character_
    )
  )

net_data <- net_data %>%
  mutate(
    Surv_1yr = factor(Surv_1yr, levels = c("No", "Yes")),
    Surv_3yr = factor(Surv_3yr, levels = c("No", "Yes")),
    Surv_5yr = factor(Surv_5yr, levels = c("No", "Yes"))
  )


net_data$Gender <- factor(net_data$Gender, levels = c(0, 1),
                          labels = c("Male", "Female"))
net_data$Smoke  <- factor(net_data$Smoke, levels = c(0, 1, 2),
                          labels = c("Never", "Ever", "Current"))
net_data$Cstage <- factor(net_data$Cstage, levels = c(1, 2, 3, 4),
                          labels = c('I', 'II', 'III', 'IV'))

net_data <- net_data %>% select(-c(OS_diag, Dead, os_year))


name_map <- setNames(data_origin$Compounds, data_origin$Index)
old_names <- colnames(net_data)
metab_idx <- which(old_names %in% names(name_map))
new_names <- old_names
new_names[metab_idx] <- name_map[old_names[metab_idx]]
new_names[which(new_names == "Kaempferol 3-O-(6''-O-acetyl)glucoside-7-O-rhamnoside")] <- 'K-3-(6"Ac)Gal-7-Rha'
new_names[which(new_names == "11beta-Hydroxyandrosterone")] <- '11βOHAn'
colnames(net_data) <- new_names

selected_metabs <- colnames(net_data)[5:13]

explist_names <- name_map[explist] 
metablist_names <- name_map[metablist]




mi_clinical_metab <- data.frame()

for (v in adjust_vars) {
  for (m in selected_metabs) {
    v_is_num <- is.numeric(net_data[[v]])
    m_is_num <- is.numeric(net_data[[m]])
    
    if (v_is_num && m_is_num) {
      test_type <- "mi-g"      
    } else {
      test_type <- "mi-cg"     
    }
    
    mi_test <- ci.test(v, m, data = net_data, test = test_type)
    mi_clinical_metab <- rbind(mi_clinical_metab, data.frame(
      From = v, To = m,
      mi = mi_test$statistic, p_value = mi_test$p.value,
      stringsAsFactors = FALSE))
  }
}


mi_clinical_metab$fdr <- p.adjust(mi_clinical_metab$p_value, method = "BH")
mi_clinical_metab <- mi_clinical_metab %>% arrange(p_value)

sig_clinical_metab <- mi_clinical_metab %>% filter(fdr < 0.05)


mi_clinical_surv <- data.frame()

for (v in adjust_vars) {
  for (s in surv_nodes) {
    v_is_num <- is.numeric(net_data[[v]])
    s_is_num <- is.numeric(net_data[[s]])
    
    if (v_is_num && s_is_num) {
      test_type <- "mi-g"       
    } else if (!v_is_num && !s_is_num) {
      test_type <- "mi"       
    } else {
      test_type <- "mi-cg"   
    }
    
    mi_test <- ci.test(v, s, data = net_data, test = test_type)
    mi_clinical_surv <- rbind(mi_clinical_surv, data.frame(
      From = v, To = s,
      mi = mi_test$statistic, p_value = mi_test$p.value,
      stringsAsFactors = FALSE))
  }
}



mi_clinical_surv$fdr <- p.adjust(mi_clinical_surv$p_value, method = "BH")
mi_clinical_surv <- mi_clinical_surv %>% arrange(p_value)

sig_clinical_surv <- mi_clinical_surv %>% filter(fdr < 0.05)


mi_metab_surv <- data.frame()

for (m in selected_metabs) {
  for (s in surv_nodes) {
    mi_test <- ci.test(m, s, data = net_data, test = "mi-cg")
    mi_metab_surv <- rbind(mi_metab_surv, data.frame(
      From = m, To = s,
      mi = mi_test$statistic, p_value = mi_test$p.value,
      stringsAsFactors = FALSE))
  }
}

mi_metab_surv$fdr <- p.adjust(mi_metab_surv$p_value, method = "BH")
mi_metab_surv <- mi_metab_surv %>% arrange(p_value)

sig_metab_surv <- mi_metab_surv %>% filter(fdr < 0.05)



mi_metab_surv_adj <- data.frame()

for (m in selected_metabs) {
  for (s in surv_nodes) {
    cmi_test <- ci.test(m, s, adjust_vars, data = net_data, test = "mi-cg")
    mi_metab_surv_adj <- rbind(mi_metab_surv_adj, data.frame(
      From = m, To = s,
      mi = cmi_test$statistic, p_value = cmi_test$p.value,
      stringsAsFactors = FALSE))
  }
}

mi_metab_surv_adj$fdr <- p.adjust(mi_metab_surv_adj$p_value, method = "BH")
mi_metab_surv_adj <- mi_metab_surv_adj %>% arrange(p_value)

sig_metab_surv_adj <- mi_metab_surv_adj %>% filter(fdr < 0.05)



mi_metab_surv_cstage <- data.frame()

for (m in selected_metabs) {
  for (s in surv_nodes) {
    cmi_test <- ci.test(m, s, "Cstage", data = net_data, test = "mi-cg")
    mi_metab_surv_cstage <- rbind(mi_metab_surv_cstage, data.frame(
      From = m, To = s,
      mi = cmi_test$statistic, p_value = cmi_test$p.value,
      stringsAsFactors = FALSE))
  }
}

mi_metab_surv_cstage$fdr <- p.adjust(mi_metab_surv_cstage$p_value, method = "BH")
mi_metab_surv_cstage <- mi_metab_surv_cstage %>% arrange(p_value)

sig_metab_surv_cstage <- mi_metab_surv_cstage %>% filter(fdr < 0.05)


sig_metabs_uncond  <- unique(sig_metab_surv$From)
sig_metabs_cstage  <- unique(sig_metab_surv_cstage$From)
sig_metabs_full    <- unique(sig_metab_surv_adj$From)

sig_metabs_any <- unique(c(sig_metabs_uncond, sig_metabs_cstage, sig_metabs_full))
sig_metabs_all <- Reduce(intersect, list(sig_metabs_uncond, sig_metabs_cstage, sig_metabs_full))




bl_metab <- data.frame(stringsAsFactors = FALSE)

bl_metab <- rbind(bl_metab, expand.grid(
  from = surv_nodes,
  to   = c(adjust_vars, selected_metabs),
  stringsAsFactors = FALSE))

bl_metab <- rbind(bl_metab, data.frame(
  from = c("Surv_3yr", "Surv_5yr", "Surv_5yr"),
  to   = c("Surv_1yr", "Surv_1yr", "Surv_3yr"),
  stringsAsFactors = FALSE))

bl_metab <- rbind(bl_metab, data.frame(
  from = c(selected_metabs, "Age_blood", "Smoke", "Cstage"),
  to   = rep("Gender", 12), stringsAsFactors = FALSE))

bl_metab <- rbind(bl_metab, data.frame(
  from = c(selected_metabs, "Gender", "Smoke", "Cstage"),
  to   = rep("Age_blood", 12), stringsAsFactors = FALSE))

bl_metab <- rbind(bl_metab, expand.grid(
  from = selected_metabs[5:9],
  to   = selected_metabs[1:4],
  stringsAsFactors = FALSE))

bl_metab <- rbind(bl_metab, expand.grid(
  from = selected_metabs,
  to   = c("Cstage", "Smoke"),
  stringsAsFactors = FALSE))


set.seed(1234)
boot_metab <- boot.strength(
  net_data,
  R         = 500,
  algorithm = "tabu",
  algorithm.args = list(
    score     = "bic-cg",
    blacklist = bl_metab,
    tabu      = 50,
    max.iter  = 10000))

dag_metab <- averaged.network(boot_metab, threshold = 0.5)


# ==================== constrcut net ====================

all_nodes <- nodes(dag_metab)


bootstrap_arcs <- arcs(dag_metab)
bootstrap_arcs <- bootstrap_arcs[
  !(bootstrap_arcs[, "from"] %in% surv_nodes) &
    !(bootstrap_arcs[, "to"]   %in% surv_nodes), , drop = FALSE]


dag_final <- empty.graph(all_nodes)
arcs(dag_final) <- bootstrap_arcs


if (nrow(sig_clinical_metab) > 0) {
  for (i in 1:nrow(sig_clinical_metab)) {
    dag_final <- set.arc(dag_final, sig_clinical_metab$From[i], sig_clinical_metab$To[i])
  }

}

if (nrow(sig_clinical_surv) > 0) {
  for (i in 1:nrow(sig_clinical_surv)) {
    dag_final <- set.arc(dag_final, sig_clinical_surv$From[i], sig_clinical_surv$To[i])
  }
  
}


if (nrow(sig_metab_surv_adj) > 0) {
  for (i in 1:nrow(sig_metab_surv_adj)) {
    dag_final <- set.arc(dag_final, sig_metab_surv_adj$From[i], sig_metab_surv_adj$To[i])
  }

  sig_edges_used <- sig_metab_surv_adj
  sig_label <- "Conditions based on all clinical variables"
} else if (nrow(sig_metab_surv_cstage) > 0) {
  for (i in 1:nrow(sig_metab_surv_cstage)) {
    dag_final <- set.arc(dag_final, sig_metab_surv_cstage$From[i], sig_metab_surv_cstage$To[i])
  }

  sig_edges_used <- sig_metab_surv_cstage
  sig_label <- "Conditions based on stage"
} else if (nrow(sig_metab_surv) > 0) {
  for (i in 1:nrow(sig_metab_surv)) {
    dag_final <- set.arc(dag_final, sig_metab_surv$From[i], sig_metab_surv$To[i])
  }

  sig_edges_used <- sig_metab_surv
  sig_label <- "Unconditional"
} else {
  cat("Warning: No significant metabolites → Survival branch.")
  sig_edges_used <- data.frame()
  sig_label <- "No significant edge"
}

cat("The final network has a total of", nrow(dag_final$arcs), "edges.")

# ==================== plot ====================

library(igraph)
library(ggraph)
library(ggplot2)

all_nodes_final <- nodes(dag_final)
clinical_nodes <- adjust_vars
metab_nodes <- setdiff(all_nodes_final, c(surv_nodes, clinical_nodes))

g <- igraph::graph_from_data_frame(d = dag_final$arcs, directed = TRUE)

V(g)$category <- ifelse(names(V(g)) %in% surv_nodes, "Survival",
                        ifelse(names(V(g)) %in% clinical_nodes, "Clinical",
                               ifelse(names(V(g)) %in% explist_names, "Exogenous", "Endogenous")))

edge_from <- names(V(g))[igraph::tail_of(g, igraph::E(g))]
edge_to   <- names(V(g))[igraph::head_of(g, igraph::E(g))]

E(g)$edge_type <- "Metabolite-Metabolite"
is_metab_surv <- edge_from %in% selected_metabs & edge_to %in% surv_nodes
E(g)$edge_type[is_metab_surv] <- "Metabolite→Survival"


r_inner <- 0.25
triangle_angles <- seq(90, 450, length.out = 4)[1:3]
triangle_angles_rad <- triangle_angles * pi / 180

inner_layout <- data.frame(
  x = r_inner * cos(triangle_angles_rad),
  y = r_inner * sin(triangle_angles_rad),
  name = surv_nodes,
  stringsAsFactors = FALSE
)

exogenous_in_graph  <- intersect(metab_nodes, explist_names)
endogenous_in_graph <- intersect(metab_nodes, metablist_names)

if (length(exogenous_in_graph) == 0 && length(endogenous_in_graph) == 0) {
  outer_nodes_grouped <- c(clinical_nodes, metab_nodes)
} else {
  outer_nodes_grouped <- c(clinical_nodes, exogenous_in_graph, endogenous_in_graph)
}

missing_nodes <- setdiff(c(clinical_nodes, metab_nodes), outer_nodes_grouped)
if (length(missing_nodes) > 0) {
  outer_nodes_grouped <- c(outer_nodes_grouped, missing_nodes)
}

n_outer <- length(outer_nodes_grouped)

outer_angles <- seq(0, 360, length.out = n_outer + 1)[1:n_outer]
outer_angles_rad <- outer_angles * pi / 180

outer_layout <- data.frame(
  x = 1.0 * cos(outer_angles_rad),
  y = 1.0 * sin(outer_angles_rad),
  name = outer_nodes_grouped,
  stringsAsFactors = FALSE
)

layout_df <- rbind(inner_layout, outer_layout)
layout_df <- layout_df[match(names(V(g)), layout_df$name), ]
layout_matrix <- as.matrix(layout_df[, c("x", "y")])

node_colors <- c(
  "Survival"   = "#EBA65D",
  "Clinical"   = "#B7D07E",
  "Exogenous"  = "#91D9E7",   
  "Endogenous" = "#8FA5D5"    
)

edge_colors <- c(
  "Metabolite→Survival"    = "#CA3730",
  "Metabolite-Metabolite"  = "grey60"
)

edge_alphas <- c(
  "Metabolite→Survival"    = 0.8,
  "Metabolite-Metabolite"  = 0.6
)

node_angle <- function(x, y) {
  angle <- atan2(y, x) * 180 / pi
  angle <- ifelse(angle < 0, angle + 360, angle)
  return(angle)
}

node_hjust <- function(x, y) {
  angle <- atan2(y, x) * 180 / pi
  angle <- ifelse(angle < 0, angle + 360, angle)
  hjust <- ifelse(angle > 90 & angle < 270, 1.4, -0.4)
  return(hjust)
}

set.seed(1105)

p <- ggraph(g, layout = "manual", x = layout_matrix[, 1], y = layout_matrix[, 2]) +
  geom_edge_arc(
    aes(color = edge_type, alpha = edge_type,  edge_width = edge_type),
    arrow = arrow(length = unit(1.5, 'mm'), type = "closed", angle = 20),
    end_cap = circle(5, 'mm'),
    start_cap = circle(5, 'mm'),
    strength = 0.15
  ) +
  geom_node_point(
    aes(color = category),
    size = 8,
    fill = "white",
    stroke = 1.8
  ) +
  geom_node_text(
    aes(
      label = name,
      angle = 0,
      hjust = node_hjust(x, y)
    ),
    size = 3.2,
    fontface = "plain"
  ) +
  scale_edge_width_manual(
    values = c(
      "Metabolite→Survival"    = 0.5,
      "Metabolite-Metabolite"  = 0.4
    ),
    guide = "none"
  )+
  scale_color_manual(values = node_colors, guide = "none") +
  scale_edge_color_manual(
    values = c("Metabolite→Survival" = "#CA3730", "Metabolite-Metabolite" = "grey60"),
    guide = "none"
  ) +
  scale_edge_alpha_manual(values = edge_alphas, guide = "none") +
  # 图例
  annotate("point", x = -0.8, y = -1.205, color = "#EBA65D", size = 4) +
  annotate("text",  x = -0.73, y = -1.2, label = "Survival", hjust = 0, size = 3.5) +
  annotate("point", x = -0.35, y = -1.205, color = "#91D9E7", size = 4) +
  annotate("text",  x = -0.27, y = -1.2, label = "Exogenous", hjust = 0, size = 3.5) +
  annotate("point", x = 0.2, y = -1.205, color = "#8FA5D5", size = 4) +
  annotate("text",  x = 0.27, y = -1.2, label = "Endogenous", hjust = 0, size = 3.5) +
  annotate("point", x = 0.78, y = -1.205, color = "#B7D07E", size = 4) +
  annotate("text",  x = 0.85, y = -1.2, label = "Clinical", hjust = 0, size = 3.5) +
  coord_fixed(clip = "off") +
  theme_graph(base_size = 14, base_family = "Arial") +
  theme(plot.margin = margin(10, 50, 10, 50))

print(p)

ggsave("03result/01Figure/Fig3B.pdf", p, width = 6, height = 6, device = cairo_pdf)


