rm(list = ls())
# setwd('D:/Cancer Model/03LC/01LC_Metab/v3/')

library(dplyr)
library(tidyr)
library(survival)
library(survminer)
library(timeROC)

load('LCmetab-main/data/demo_BLCS_data.RData')

metab_compute <- demo_BLCS_data
################AUC#######################################
betalist_A <- summary(coxph(Surv(OS_diag,Dead)~Smoke+Ctype+Cstage+Gender+Age_blood+
                              MW0138609+MW0169673+
                              MW0105314+MW0049035+MW0012126+MW0016406+
                              MW0011108*Age_blood+MW0062208*Age_blood+
                              MEDP1685*Age_blood,
                            data=metab_compute))$coefficients[,1]

betalist_C <- summary(coxph(Surv(OS_diag,Dead)~Smoke+Ctype+Cstage+Gender+Age_blood,
                            data=metab_compute))$coefficients[1]

model_A <- coxph(Surv(OS_diag,Dead)~Smoke+Ctype+Cstage+Gender+Age_blood+
                   MW0138609+MW0169673+
                   MW0105314+MW0049035+MW0012126+MW0016406+
                   MW0011108*Age_blood+MW0062208*Age_blood+
                   MEDP1685*Age_blood,
                          data=metab_compute)
model_C <- coxph(Surv(OS_diag,Dead)~Smoke+Ctype+Cstage+Gender+Age_blood,
                 data=metab_compute)

metab_compute$All_covar <- predict(model_A,type="lp",newdata = metab_compute)
metab_compute$Clin_covar <- predict(model_C,type="lp",newdata = metab_compute)

ROC_A <- timeROC(T=metab_compute$OS_diag,
                 delta=metab_compute$Dead,marker=metab_compute$All_covar,
                 cause=1,weighting="marginal",iid=T,
                 times=c(1,2,3,4,5,6,7,8,9,10))

ROC_C <- timeROC(T=metab_compute$OS_diag,
                 delta=metab_compute$Dead,marker=metab_compute$Clin_covar,
                 cause=1,weighting="marginal",iid=T,
                 times=c(1,2,3,4,5,6,7,8,9,10))

AUC_matrix <- data.frame(AUC_A=ROC_A$AUC,
                         AUC_C=ROC_C$AUC,
                         CIL_A=confint(ROC_A)$CI_AUC[,1]/100,
                         CIU_A=confint(ROC_A)$CI_AUC[,2]/100,
                         CIL_C=confint(ROC_C)$CI_AUC[,1]/100,
                         CIU_C=confint(ROC_C)$CI_AUC[,2]/100,
                         P=compare(ROC_A,ROC_C)$p_values_AUC)

AUC_matrix$P <- format(AUC_matrix$P, scientific = TRUE, digits = 3)

AUC_matrix$time <- c(1:10)


library(ggplot2)
library(tidyr)

data_long <- pivot_longer(
  AUC_matrix,
  cols = starts_with(c("AUC_", "CIL_", "CIU_")), 
  names_to = c(".value", "Model_ID"),             
  names_sep = "_"                                
)
data_long$Model <- ifelse(data_long$Model_ID == 'A', 'Metabolite-augmented', 'Baseline')
data_long <- select(data_long, c('time','Model','Model_ID','AUC','CIL','CIU','P'))


bar_width <- 0.7
dodge_width <- 0.7
AUC_matrix <- AUC_matrix %>%
  mutate(
    x_start = time - dodge_width / 4, # Position for 'Clinical+Metabolite'
    x_end = time + dodge_width / 4,   # Position for 'Clinical'
    y_pos = pmax(AUC_A, AUC_C)+0.02 , # Position horizontal line above the taller bar
    label_y = y_pos + 0.005 # Position for the significance label
  )

# Plotting with ggplot2
p1 <- ggplot(data_long, aes(x = factor(time), y = AUC, fill = Model)) + 
  geom_col(position = position_dodge(width = dodge_width), width = bar_width) +
  # Add the significance brackets and labels using the prepared AUC_matrix
  geom_segment(data = AUC_matrix, aes(x = x_start, xend = x_end, y = max(AUC_A)+0.003, yend = max(AUC_A)+0.003), linewidth = 0.3, inherit.aes = FALSE) +
  geom_segment(data = AUC_matrix, aes(x = x_start, xend = x_start, y = max(AUC_A)+0.003, yend = AUC_C+0.001), linewidth = 0.3, inherit.aes = FALSE) +
  geom_segment(data = AUC_matrix, aes(x = x_end, xend = x_end, y = AUC_A+0.001, yend = max(AUC_A)+0.003), linewidth = 0.3, inherit.aes = FALSE) +
  geom_text(
    data = AUC_matrix,
    aes(x = factor(time), y = max(AUC_A)+0.011, label = P),
    inherit.aes = FALSE,
    size = 3
  ) +
  scale_fill_manual(values = c("Metabolite-augmented" = "#FBAFA0", "Baseline" = "#6D9FC6")) + 
  labs( 
    title = 'Discrimination ability', 
    x = "Time (years)", 
    y = "Time-dependent AUC", 
    fill = "Model" 
  ) + 
  scale_y_continuous(expand = c(0.0,0), breaks = seq(from = 0.7, to = 0.9, by = 0.05))+
  coord_cartesian(ylim = c(0.74, 0.89)) + # Increased upper limit and turned clipping off
  theme_bw() + 
  theme( 
    panel.grid.minor.y = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"), 
    legend.position = 'bottom', 
    plot.caption = element_text(hjust = 1, face = "italic"), 
    text = element_text(size = 12),
    # Add plot.margin to ensure there's space for labels on top
    plot.margin = margin(t = 20, r = 5.5, b = 5.5, l = 5.5, unit = "pt")
  )

pdf(file = '03result/01Figure/Fig4A_AUC.pdf',width = 8,height = 5)
p1
dev.off()

################ Optimistic Bias ###############################
betalist_A <- summary(coxph(Surv(OS_diag,Dead)~Smoke+Ctype+Cstage+Gender+Age_blood+
                              MW0138609+MW0169673+
                              MW0105314+MW0049035+MW0012126+MW0016406+
                              MW0011108*Age_blood+MW0062208*Age_blood+
                              MEDP1685*Age_blood,
                            data=metab_compute))$coefficients[,1]

betalist_C <- summary(coxph(Surv(OS_diag,Dead)~Smoke+Ctype+Cstage+Gender+Age_blood,
                            data=metab_compute))$coefficients[1]

# ----------------- 10-fold Cross-Validation -----------------
set.seed(1091) 
k_folds <- 10


fold_ids <- sample(rep(1:k_folds, length.out = nrow(metab_compute)))

cv_All_covar <- numeric(nrow(metab_compute))
cv_Clin_covar <- numeric(nrow(metab_compute))

for (i in 1:k_folds) {
  train_idx <- which(fold_ids != i) # train set (K-1)
  test_idx <- which(fold_ids == i)  # validation set (1)
  
  fit_A_cv <- coxph(Surv(OS_diag,Dead)~Smoke+Ctype+Cstage+Gender+Age_blood+
                      MW0138609+MW0169673+MW0105314+MW0049035+MW0012126+MW0016406+
                      MW0011108*Age_blood+MW0062208*Age_blood+MEDP1685*Age_blood,
                    data = metab_compute[train_idx, ])
  
  fit_C_cv <- coxph(Surv(OS_diag,Dead)~Smoke+Ctype+Cstage+Gender+Age_blood,
                    data = metab_compute[train_idx, ])
  
  cv_All_covar[test_idx] <- predict(fit_A_cv, newdata = metab_compute[test_idx, ], type = "lp")
  cv_Clin_covar[test_idx] <- predict(fit_C_cv, newdata = metab_compute[test_idx, ], type = "lp")
}


metab_compute$All_covar <- cv_All_covar
metab_compute$Clin_covar <- cv_Clin_covar

# ----------------- Cross-Validation timeROC -----------------
ROC_A <- timeROC(T=metab_compute$OS_diag,
                 delta=metab_compute$Dead,marker=metab_compute$All_covar,
                 cause=1,weighting="marginal",iid=T,
                 times=c(1,2,3,4,5,6,7,8,9,10))

ROC_C <- timeROC(T=metab_compute$OS_diag,
                 delta=metab_compute$Dead,marker=metab_compute$Clin_covar,
                 cause=1,weighting="marginal",iid=T,
                 times=c(1,2,3,4,5,6,7,8,9,10))

AUC_matrix <- data.frame(AUC_A=ROC_A$AUC,
                         AUC_C=ROC_C$AUC,
                         CIL_A=confint(ROC_A)$CI_AUC[,1]/100,
                         CIU_A=confint(ROC_A)$CI_AUC[,2]/100,
                         CIL_C=confint(ROC_C)$CI_AUC[,1]/100,
                         CIU_C=confint(ROC_C)$CI_AUC[,2]/100,
                         P=compare(ROC_A,ROC_C)$p_values_AUC)

AUC_matrix$P <- format(AUC_matrix$P, scientific = TRUE, digits = 3)

AUC_matrix$time <- c(1:10)

library(ggplot2)
library(tidyr)

data_long <- pivot_longer(
  AUC_matrix,
  cols = starts_with(c("AUC_", "CIL_", "CIU_")), 
  names_to = c(".value", "Model_ID"),             
  names_sep = "_"                                
)
data_long$Model <- ifelse(data_long$Model_ID == 'A', 'Metabolite-augmented', 'Baseline')
data_long <- select(data_long, c('time','Model','Model_ID','AUC','CIL','CIU','P'))

write.csv(data_long, '03result/02Table/TabS13_AUC_compare_optbia.csv',row.names = F)


################E/O ratio#######################################
library(metamisc)
model_A <- coxph(Surv(OS_diag,Dead)~All_covar,
                 data=metab_compute)
model_C <- coxph(Surv(OS_diag,Dead)~Clin_covar,
                 data=metab_compute)

fit1 <- survfit(Surv(OS_diag,Dead)~1,metab_compute)
fit2 <- survfit(coxph(Surv(OS_diag,Dead) ~ Clin_covar,data = metab_compute))
fit3 <- survfit(coxph(Surv(OS_diag,Dead) ~ All_covar,data = metab_compute))

ds <- metab_compute
plotname_list <- c(2:3)
times <- c(1:10)

res <- data.frame()
for (i in 1:length(plotname_list)) {
  Model <- ifelse(i==2,'Baseline','Metabolite-augmented')
  fit1 <- survfit(Surv(OS_diag,Dead)~1,ds)
  model <- get(paste0('fit',plotname_list[i]))
  
  for (j in 1:length(times)) {
    time <- times[j]
    Po <- summary(fit1,time= time)$surv
    Pe <- summary(model,time = time)$surv
    est_t = oecalc(Pe = Pe,  Po = Po, N = nrow(ds))
    
    EO_ratio_t = 1/est_t$theta
    
    EO_ratio_LCI_t = 1/est_t$theta.ciub
    EO_ratio_UCI_t = 1/est_t$theta.cilb
    resi <- data.frame(Model=Model,time=time,
                       EO_ratio=EO_ratio_t,EO_ratio_LCI=EO_ratio_LCI_t,EO_ratio_UCI=EO_ratio_UCI_t)
    res <- rbind(res,resi)
  }
}
EOlist <- res

EOlist <- select(EOlist, c('time','Model','EO_ratio',"EO_ratio_LCI","EO_ratio_UCI")) %>%
  arrange(time)
write.csv(EOlist, '03result/02Table/EO_compare.csv',row.names = F)

p2 <- ggplot(EOlist, aes(x = factor(time), y = EO_ratio, color = Model)) + # 1. 交換 x 和 y
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  
  geom_errorbar(aes(ymin = EO_ratio_LCI, ymax = EO_ratio_UCI), 
                position = position_dodge(width = 0.7), width = 0.2) +
  
  geom_hline(yintercept = 1, linetype = "dashed") +
  
  scale_y_continuous(name = "E/O ratio") +
  scale_color_manual(values = c("Metabolite-augmented" = "#a52a2a", 'Baseline'="#1f78b4")) +
  labs(color = "Model", 
       x = "Time (years)", 
       y = "E/O ratio",
       title = "Calibration performance") + 
  theme_bw() +
  theme(
    legend.position = "bottom", 
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )

print(p2)

pdf(file = '03result/01Figure/Fig4B_EO.pdf',width = 8,height = 5)
p2
dev.off()


#################### DCA ######################################################

source('D:/Cancer Model/my_function/stdca.R')

options(scipen = 0)
Alldata <- na.omit(metab_compute)
model_C <- coxph(Surv(OS_diag,Dead)~Clin_covar,data = Alldata)
model_A <- coxph(Surv(OS_diag,Dead)~All_covar,data = Alldata)

Alldata <- Alldata[complete.cases(Alldata[,c('Clin_covar','All_covar')]),]


Alldata$Covariate3 <- c(1- (summary(survfit(model_C,
                                            newdata=Alldata), times=3)$surv))
Alldata$All3 <- c(1- (summary(survfit(model_A,
                                          newdata=Alldata), times=3)$surv))

Alldata$Covariate5 <- c(1- (summary(survfit(model_C,
                                            newdata=Alldata), times=5)$surv))
Alldata$All5 <- c(1- (summary(survfit(model_A,
                                          newdata=Alldata), times=5)$surv))


pdf(file="03result/01Figure/Fig4C_DCA.pdf", width=7.5, height=7.5)
par(mfrow=c(2,2), mgp=c(2.5,1,0), mar=c(4,4,2,0.5))

xstart <- 0.2
threshold <- 0.3

ymin <- -0.1
fit1 <- stdca(data=Alldata, outcome="Dead", ttoutcome="OS_diag", timepoint=3,
              predictors=c('Covariate3','All3'), xstart=xstart, xstop=0.8,
              intervention=FALSE, xby=0.01, ymin=ymin)
idx <- which.min(abs(fit1$net.benefit$threshold - threshold))
y_vals <- sapply(fit1$predictors$predictor, function(p) fit1$net.benefit[idx, p])
max_y <- max(y_vals, na.rm=TRUE)
segments(x0=threshold, y0=ymin, x1=threshold, y1=max_y, lty=2, col="darkgray")
for (pred in fit1$predictors$predictor) {
  y <- fit1$net.benefit[idx, pred]
  segments(x0=xstart, y0=y, x1=threshold, y1=y, lty=2, col="darkgray")
  text(x=xstart, y=y+0.01, labels=round(y, 3), cex=1, adj=0, col="black")
}


ymin <- 0
fit2 <- stdca(data=Alldata, outcome="Dead", ttoutcome="OS_diag", timepoint=3,
              predictors=c('Covariate3','All3'), xstart=xstart, xstop=0.8,
              intervention=TRUE, xby=0.01)
idx <- which.min(abs(fit2$interventions.avoided$threshold - threshold))
y_vals <- sapply(fit2$predictors$predictor, function(p) fit2$interventions.avoided[idx, p])
max_y <- max(y_vals, na.rm=TRUE)
segments(x0=threshold, y0=ymin, x1=threshold, y1=max_y, lty=2, col="darkgray")
for (pred in fit2$predictors$predictor) {
  y <- fit2$interventions.avoided[idx, pred]
  segments(x0=xstart, y0=y, x1=threshold, y1=y, lty=2, col="darkgray")
  text(x=xstart, y=y+0.01, labels=round(y, 3), cex=1, adj=0, col="black")
}

ymin <- -0.1
fit3 <- stdca(data=Alldata, outcome="Dead", ttoutcome="OS_diag", timepoint=5,
              predictors=c('Covariate5','All5'), xstart=xstart, xstop=0.8,
              intervention=FALSE, xby=0.01, ymin=ymin)
idx <- which.min(abs(fit3$net.benefit$threshold - threshold))
y_vals <- sapply(fit3$predictors$predictor, function(p) fit3$net.benefit[idx, p])
max_y <- max(y_vals, na.rm=TRUE)
segments(x0=threshold, y0=ymin, x1=threshold, y1=max_y, lty=2, col="darkgray")
for (pred in fit3$predictors$predictor) {
  y <- fit3$net.benefit[idx, pred]
  segments(x0=xstart, y0=y, x1=threshold, y1=y, lty=2, col="darkgray")
  text(x=xstart, y=y+0.01, labels=round(y, 3), cex=1, adj=0, col="black")
}

ymin <- 0
fit4 <- stdca(data=Alldata, outcome="Dead", ttoutcome="OS_diag", timepoint=5,
              predictors=c('Covariate5','All5'), xstart=xstart, xstop=0.8,
              intervention=TRUE, xby=0.01)
idx <- which.min(abs(fit4$interventions.avoided$threshold - threshold))
y_vals <- sapply(fit4$predictors$predictor, function(p) fit4$interventions.avoided[idx, p])
max_y <- max(y_vals, na.rm=TRUE)
segments(x0=threshold, y0=ymin, x1=threshold, y1=max_y, lty=2, col="darkgray")
for (pred in fit4$predictors$predictor) {
  y <- fit4$interventions.avoided[idx, pred]
  segments(x0=xstart, y0=y, x1=threshold, y1=y, lty=2, col="darkgray")
  text(x=xstart, y=y+0.01, labels=round(y, 3), cex=1, adj=0, col="darkgray")
}

dev.off()


# 3-year
a <- fit1$net.benefit
# a <- a[a$threshold >=0 & a$threshold <=0.5,]
a$delta <- a$All3- a$Covariate3
mean(a$All3);mean(a$Covariate3)
t.test(a$All3,a$Covariate3,paired = T)$p.value
mean(a$delta)

b <- fit1$interventions.avoided
# b <- b[b$threshold >=0 & b$threshold <=0.5,]
b$delta <- b$All3- b$Covariate3
mean(b$All3);mean(b$Covariate3)
t.test(b$All3,b$Covariate3,paired = T)$p.value
mean(b$delta)

# 5-year
a <- fit3$net.benefit
# a <- a[a$threshold >=0 & a$threshold <=0.5,]
a$delta <- a$All5- a$Covariate5
mean(a$All5);mean(a$Covariate5)
t.test(a$All5,a$Covariate5,paired = T)$p.value
mean(a$delta)

b <- fit3$interventions.avoided
# b <- b[b$threshold >=0 & b$threshold <=0.5,]
b$delta <- b$All5- b$Covariate5
mean(b$All5);mean(b$Covariate5)
t.test(b$All5,b$Covariate5,paired = T)$p.value
mean(b$delta)


#############NRI_IDI##################################
library(survIDINRI)

z.std <- as.matrix(subset(metab_compute, select = c('Clin_covar')))  
z.new <- as.matrix(subset(metab_compute, select = c('All_covar')))   

times <- 1:10

results_list <- list()

for (t in times) {
  cat("Processing time =", t, "years\n")
  
  # caculate IDI/NRI
  fit <- IDI.INF(
    indata = metab_compute[, c('OS_diag', 'Dead')],  
    covs0 = z.std,
    covs1 = z.new,
    t0 = t,
    npert = 1000
  )
  
  out <- IDI.INF.OUT(fit) 
  df <- as.data.frame(out)
  df$Time <- t
  df$Metric <- c('IDI','NRI','Median improvement')  # M1 = IDI, M2 = NRI, M3 = ADI
  results_list[[paste0("t", t)]] <- df
}


all_results <- do.call(rbind, results_list)
rownames(all_results) <- NULL


all_results <- all_results[, c("Time", "Metric", "Est.", "Lower", "Upper", "p-value")]


write.csv(all_results, file = "03result/02Table/TabS14_IDI_NRI_1to10yr.csv", row.names = FALSE)

# 查看前几行
print(all_results)
