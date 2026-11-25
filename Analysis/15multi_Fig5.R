rm(list=ls())
gc() 

library(tidyverse);library(data.table)

library(readxl)
library(survival)
metab_protein <- read_xlsx('/data1/user/mjxue/05UKB_LC/metab2protein.xlsx',sheet = 'metab2protein')

load('/data1/user/mjxue/05UKB_LC/ukb_protein_analysis.RData')

ukb_genome <- fread("/data1/user/anonymous/Olink/LC_HSPH_Olink_Proteomics_analysis_20250304/data/UKB_genome_predicted.csv")
gene_anno <- fread('/data1/user/xxwen/eQTL_Gen/annotation/annotation.txt')

ukb_genome <- data.frame(ukb_genome)

gene_sel <- paste0(gene_anno$SiteName[gene_anno$GeneChr %in% c(1:22)],'.ma')

ukb_genome <- ukb_genome[,which(colnames(ukb_genome) %in% c("family.ID","sample.ID",gene_sel))]

id_inter <- intersect(ukb_genome$sample.ID, ukb_protein_analysis$eid_chen)

ukb_genome <- ukb_genome[which(ukb_genome$sample.ID %in% id_inter),]
ukb_pro <- ukb_protein_analysis[which(ukb_protein_analysis$eid_chen %in% id_inter),]

merged_data <- left_join(ukb_genome,ukb_pro,by = c("family.ID"="eid_chen"))


################ protein-prognosis ###################################
inter_pro <- intersect(colnames(ukb_pro),metab_protein$Gene)

res <- data.frame()
for (i in 1:length(inter_pro)) {
  protein <- inter_pro[i]
  fit <- coxph(as.formula(paste0('Surv(OS_time, Dead) ~ ', protein,' + age+sex+as.factor(smokeingstatus)+as.factor(histo_new)')),
               data = ukb_protein_analysis)
  
  HR <- summary(fit)$coefficients[protein, 'exp(coef)']
  LCI <- summary(fit)$conf.int[protein, 'lower .95']
  UCI <- summary(fit)$conf.int[protein, 'upper .95']
  P <- summary(fit)$coefficients[protein, 'Pr(>|z|)']
  resi <- data.frame(protein, HR, LCI, UCI, P)
  res <- rbind(res,resi)
}

res[res$P<=0.05,]

pro_sur_res <- res
write.csv(pro_sur_res, '05UKB_LC/multiomic/TabS16pro_sur_res.csv', row.names = F)

################ gene-cor-protein ###################################
library(foreach)
library(doParallel)
library(survival)

num_cores <- 15
cl <- makeCluster(num_cores)
registerDoParallel(cl)


genelist <- colnames(ukb_genome)[-1:-2]
res_gene_parallel <- foreach(
  gene = genelist, 
  .combine = 'rbind', 
  .packages = 'survival'
) %dopar% {
  
  
  formula_str <- paste0('Surv(OS_time, Dead) ~ ', gene, ' + age+sex+as.factor(smokeingstatus)+as.factor(histo_new)')
  fit <- coxph(as.formula(formula_str), data = merged_data)
  
  HR <- summary(fit)$coefficients[1, 'exp(coef)']
  CIL <- summary(fit)$conf.int[1, 'lower .95']
  CIU <- summary(fit)$conf.int[1, 'upper .95']
  P <- summary(fit)$coefficients[1, 'Pr(>|z|)']
  
  
  data.frame(gene = gene, HR = HR, CIL = CIL, CIU = CIU, P = P)
}


stopCluster(cl)

write.csv(res_gene_parallel,'multiomic/res_gene_cox.csv',row.names = F)


which(colnames(merged_data) == "CNDP1") # 13249


# partial cor ####
library(parallel)
library(data.table)

gene_cols    <- 3:12635
protein_cols <- c(13249)

cl <- makeCluster(10)
clusterExport(cl, c("merged_data","gene_cols","protein_cols"))
clusterEvalQ(cl, library(stats))  

res_list <- parLapply(cl, gene_cols, function(gcol) {
  gene <- colnames(merged_data)[gcol]
  
  do.call(rbind, lapply(protein_cols, function(pcol){
    prot <- colnames(merged_data)[pcol]
    
    fit1 <- lm(
      as.formula(paste0(
        prot, " ~ ", gene,
        " + age+sex+as.factor(smokeingstatus)+as.factor(histo_new)"
      )),
      data = merged_data
    )
    
    fit_gene <- lm(
      as.formula(paste0(
        gene, " ~ age+sex+as.factor(smokeingstatus)+as.factor(histo_new)"
      )),
      data = merged_data
    )
    fit_prot <- lm(
      as.formula(paste0(
        prot, " ~ age+sex+as.factor(smokeingstatus)+as.factor(histo_new)"
      )),
      data = merged_data
    )
    
    r_raw <- cor(merged_data[[gene]], merged_data[[prot]],
                 use = "pairwise.complete.obs")
    pcor  <- cor(residuals(fit_gene), residuals(fit_prot))
    
    data.frame(
      gene    = gene,
      protein = prot,
      beta    = coef(summary(fit1))[gene, "Estimate"],
      se      = coef(summary(fit1))[gene, "Std. Error"],
      pval    = coef(summary(fit1))[gene, "Pr(>|t|)"],
      r_raw   = r_raw,
      pcor    = pcor,
      stringsAsFactors = FALSE
    )
  }))
})

stopCluster(cl)

final_df <- do.call(rbind, res_list)
write.csv(final_df,file = "multiomic/gene_prot_cor.csv",row.names = F)

library(data.table)
res_gene_cox <- fread('multiomic/res_gene_cox.csv')
cor_trans_pro <- read.csv("multiomic/gene_prot_cor.csv")

cor_trans_pro <- cor_trans_pro[which(cor_trans_pro$gene %in% res_gene_cox$gene[res_gene_cox$P<=0.05]),]


cor_trans_pro <- cor_trans_pro %>%
  group_by(protein) %>%
  mutate(
    p_fdr = p.adjust(pval, method = "BH")
  ) %>%
  ungroup()


cor_trans_pro[cor_trans_pro$p_fdr<=0.05,]

################ methy-gene ###################################
rm(list=ls())
gc() 

library(dplyr)
library(data.table)

# data load ####
load('/home/Users/mjxue/LC/ukb_protein_analysis.RData')
ukb_epigenome <- fread("/home/Users/mjxue/LC/UKB_epigenome_predicted.csv")
ukb_genome <- fread("/home/Users/mjxue/LC/UKB_genome_predicted.csv")
epic_anno <- fread('/home/Users/mjxue/LC/humanmethylation450_15017482_v1-2_anno.csv')


ukb_epigenome <- data.frame(ukb_epigenome)

CpG_sel <- paste0(epic_anno$IlmnID[epic_anno$CHR %in% c(1:22)],'.ma')

ukb_epigenome <- ukb_epigenome[,which(colnames(ukb_epigenome) %in% c("family.ID","sample.ID",CpG_sel))]

id_inter <- intersect(ukb_genome$sample.ID, ukb_protein_analysis$eid_chen)

ukb_epigenome <- ukb_epigenome[which(ukb_epigenome$sample.ID %in% id_inter),]
ukb_genome <- ukb_genome[which(ukb_genome$sample.ID %in% id_inter),]
ukb_pro <- ukb_protein_analysis[which(ukb_protein_analysis$eid_chen %in% id_inter),]

merged_data <- merge(ukb_epigenome,ukb_genome)
merged_data <- left_join(merged_data,ukb_pro,by = c("family.ID"="eid_chen"))


# partial cor ####
library(parallel)

which(colnames(merged_data) == "ENSG00000133612.ma") # 132865

gene_cols <- 132865

epigenome_cols <- 3:129285

cl <- makeCluster(50)
clusterExport(cl, c("merged_data","gene_cols","epigenome_cols"))
clusterEvalQ(cl, library(stats)) 

res_list <- parLapply(cl, epigenome_cols, function(ecol) {
  CpG <- colnames(merged_data)[ecol]
  
  do.call(rbind, lapply(gene_cols, function(tcol){
    gene <- colnames(merged_data)[tcol]
    
    
    fit1 <- lm(
      as.formula(paste0(
        gene, " ~ ", CpG,
        " + age+sex+as.factor(smokeingstatus)+as.factor(histo_new)"
      )),
      data = merged_data
    )
    
    
    fit_CpG <- lm(
      as.formula(paste0(
        CpG, " ~ age+sex+as.factor(smokeingstatus)+as.factor(histo_new)"
      )),
      data = merged_data
    )
    fit_gene <- lm(
      as.formula(paste0(
        gene, " ~ age+sex+as.factor(smokeingstatus)+as.factor(histo_new)"
      )),
      data = merged_data
    )
    
    
    r_raw <- cor(merged_data[[CpG]], merged_data[[gene]],
                 use = "pairwise.complete.obs")
    
    pcor  <- cor(residuals(fit_CpG), residuals(fit_gene))
    
    coef_summary <- coef(summary(fit1))
    if(CpG %in% rownames(coef_summary)) {
      beta <- coef_summary[CpG, "Estimate"]
      se <- coef_summary[CpG, "Std. Error"]
      pval <- coef_summary[CpG, "Pr(>|t|)"]
    } else {
      beta <- NA
      se <- NA
      pval <- NA
    }
    
    
    data.frame(
      CpG     = CpG,
      gene    = gene,
      beta    = beta,
      se      = se,
      pval    = pval,
      r_raw   = r_raw,
      pcor    = pcor,
      stringsAsFactors = FALSE
    )
  }))
})


stopCluster(cl)


res_df <- do.call(rbind, res_list)
write.csv(res_df,file = "/home/Users/mjxue/LC/epic_gen_cor.csv",row.names = F)

library(data.table)
res_epic_cox <- fread('/home/Users/mjxue/LC/res_epi_cox.csv')
cor_epic_trans <- read.csv("/home/Users/mjxue/LC/epic_gen_cor.csv")

cor_epic_trans <- cor_epic_trans[which(cor_epic_trans$CpG %in% res_epic_cox$gene[res_epic_cox$P<=0.05]),]



cor_epic_trans <- cor_epic_trans %>%
  group_by(gene) %>%
  mutate(
    p_fdr = p.adjust(pval, method = "BH")
  ) %>%
  ungroup()


cor_epic_trans[cor_epic_trans$p_fdr<=0.05,]


############################### plot figure ##########################################################
rm(list = ls())
library(readxl)
metab_protein <- read_xlsx('~/Cancer/LC/Metabolomics/protein_metab.xlsx',sheet = 'protein_metab')
res_pro_sur <- fread('multiomic/pro_sur_res.csv')

res_trans_cox <- fread('multiomic/res_gene_cox.csv')
cor_trans_pro <- read.csv("multiomic/gene_prot_cor.csv")

res_epic_cox <- fread('multiomic/res_epi_cox.csv')
cor_epic_trans <- read.csv("multiomic/epic_gen_cor.csv")

gene_anno <- fread('/data1/user/xxwen/eQTL_Gen/annotation/annotation.txt')
epic_anno <- fread('/data1/user/mjxue/UKB/multiomic/humanmethylation450_15017482_v1-2_anno.csv')

colnames(metab_protein)[8] <- 'protein'

cor_trans_pro <- cor_trans_pro[cor_trans_pro$gene %in% res_trans_cox$gene[res_trans_cox$P<=0.05],]
cor_trans_pro <- cor_trans_pro %>%
  group_by(protein) %>%
  mutate(
    p_fdr = p.adjust(pval, method = "BH")
  ) %>%
  ungroup()

plyr::count(cor_trans_pro$p_fdr<=0.05)
# x freq
# 1 FALSE 1441
# 2  TRUE    1

cor_epic_trans <- cor_epic_trans[cor_epic_trans$CpG %in% res_epic_cox$gene[res_epic_cox$P<=0.05],]
cor_epic_trans <- cor_epic_trans %>%
  group_by(gene) %>%
  mutate(
    p_fdr = p.adjust(pval, method = "BH")
  ) %>%
  ungroup()

plyr::count(cor_epic_trans$p_fdr<=0.05)
# x freq
# 1 FALSE 7531
# 2  TRUE   27

cor_pro_trans_sel <- cor_trans_pro[cor_trans_pro$p_fdr<=0.05,c('gene','protein','r_raw','pcor')]
cor_trans_epic_sel <- cor_epic_trans[cor_epic_trans$p_fdr<=0.05,c('CpG','gene','r_raw','pcor')]
colnames(cor_pro_trans_sel)[-1:-2] <- paste0(colnames(cor_pro_trans_sel)[-1:-2],'_progene')
colnames(cor_trans_epic_sel)[-1:-2] <- paste0(colnames(cor_trans_epic_sel)[-1:-2],'_geneepic')

ukb_pro_inter <- metab_protein[metab_protein$protein %in% res_pro_sur$protein,]

res_pro_sur_sel <- res_pro_sur[res_pro_sur$P<=0.05,]

cor_metab_pro <- left_join(res_pro_sur_sel,subset(metab_protein, select = c('Compounds','Origin','protein','Type')))

cor_metab_pro_trans <- left_join(cor_metab_pro, cor_pro_trans_sel)

cor_metab_pro_trans_epic <- full_join(cor_metab_pro_trans, cor_trans_epic_sel)

cor_metab_pro_trans_epic <- subset(cor_metab_pro_trans_epic, select = c('Compounds','protein','gene','CpG',
                                                                        'Type','HR',
                                                                        'pcor_progene',
                                                                        'pcor_geneepic'))

cor_metab_pro_trans_epic$gene <- gsub('.ma','',cor_metab_pro_trans_epic$gene)
cor_metab_pro_trans_epic$CpG <- gsub('.ma','',cor_metab_pro_trans_epic$CpG)


library(ggalluvial)
library(ggplot2)
library(dplyr)

alluvial_data <- cor_metab_pro_trans_epic %>%
  select(Compounds, protein, gene, CpG, pcor_progene, pcor_geneepic) %>%
  
  mutate(id = row_number()) %>%
  
  pivot_longer(cols = c(Compounds, protein, gene, CpG),
               names_to = "layer",
               values_to = "node") %>%
  
  filter(!is.na(node)) %>%
  
  mutate(
    layer = factor(layer, 
                   levels = c("Compounds", "protein", "gene", "CpG"),
                   labels = c("Metabolite", "Protein", "Gene", "Epigenetic"))
  ) %>%
  
  mutate(
    correlation = case_when(
      layer == "Protein" ~ pcor_progene,
      layer == "Gene" ~ pcor_geneepic,
      TRUE ~ 0
    )
  )



gg <- ggplot(alluvial_data,
             aes(x = layer, 
                 stratum = node, 
                 alluvium = id)) +
  
  geom_flow(aes(fill = correlation), 
            alpha = 0.7, 
            width = 0.3) +
  
  geom_stratum(fill = "white", 
               # alpha = 0.8, 
               color = 'black',
               width = 0.5) +
  
  geom_text(stat = "stratum", 
            aes(label = node), 
            size = 2.5) +
  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0,
                       limits = c(-1, 1),
                       name = "Partial Correlation",
                       guide = guide_colorbar(title.position = "top")) +
  scale_x_discrete(expand = c(0.05, 0.05))+
  scale_y_continuous(expand = c(0.02,0.02))+
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right",
    legend.key.width = unit(0.5, "cm")
  ) +
  labs(title = NULL)+
  annotate("segment", 
           x = -Inf, xend = Inf, 
           y = -Inf, yend = -Inf, 
           color = "black", size = 1)

print(gg)

pdf(file = 'multiomic/Fig5pathway1013.pdf', width = 8, height = 6)
gg
dev.off()
