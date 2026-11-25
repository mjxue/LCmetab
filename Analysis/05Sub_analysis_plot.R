##########################Subgroup Analysis##########################

rm(list = ls())
library(readxl)
library(dplyr)
library(survival)
setwd('D:/Cancer Model/03LC/01LC_Metab/v3/')
source('02code/function/chisq_strata_metab.R')
load('01data/02cleandata/metab_compute.RData')
load('01data/02cleandata/data_origin.RData')
load('03result/Ultra_Robust_Biomarkers.RData')

expose_list <- data_origin$Index[data_origin$Origin=='Exogenous']
metab_list <- data_origin$Index[data_origin$Origin=='Endogenous']

covar_list <- c('Age_diagnosis','Age_blood',"Gender","Smoke","Pack_year","Ctype","Cstage",'radio','chemo','bl2di','timediff')
metab_compute <- subset(metab_compute,select = c("ID",covar_list,'OS_diag','OS_bld','Dead','PFS_diag','PFS',data_origin$Index))

########################## analysis #############################
library(timeROC)
data_sub <- metab_compute
data_sub$Class_age <- ifelse(data_sub$Age_blood<=65,'Young',ifelse(is.na(data_sub$Age_blood),NA,'Older'))
data_sub$Gender <- ifelse(data_sub$Gender==0,'Female',ifelse(is.na(data_sub$Gender),NA,'Male'))
data_sub$Smoke <- ifelse(data_sub$Smoke==0,'Never',
                         ifelse(is.na(data_sub$Smoke),NA,
                                ifelse(data_sub$Smoke==2,'Former','Current')))
data_sub$Class_stage <- ifelse(data_sub$Cstage %in% c(1,2),'Early',ifelse(data_sub$Cstage %in% c(3,4),'Late',NA))
data_sub$radio <- ifelse(data_sub$radio==0,'No',ifelse(is.na(data_sub$radio),NA,'Yes'))
data_sub$chemo <- ifelse(data_sub$chemo==0,'No',ifelse(is.na(data_sub$chemo),NA,'Yes'))

coltofilter <- c("Class_age","Gender","Smoke","Ctype","Class_stage",'radio','chemo')

############ main#######################
metablist <- c(robust_exp_main$Index,robust_metab_main$Index)
HR_strat_full <- setNames(as.data.frame(matrix(nrow=1, ncol = 7)), c("Index","beta","se","pvalue","z_stat","Variable","strat_column"))
for (k in 1:length(coltofilter)) { 
  strat <- data_sub%>%pull(coltofilter[k])%>%as.factor()%>%levels
  for (i in 1:length(strat)) {
    data_sub_st <- data_sub%>%filter(!!as.symbol(coltofilter[k]) == strat[i])
    
    multi_betas <- data.frame(Index=NA,beta=NA, se=NA, pvalue=NA, z_stat=NA, number=NA, AUC3=NA,AUC5=NA,AUC10=NA)
    
    for (m in 1:length(metablist)){
      multi_betas[m,"Index"] <- metablist[m]
      metab <- metablist[m]
      #print(i)
      if (coltofilter[k]=="Class_age") {
        model <- coxph(as.formula(paste0("Surv(OS_diag,Dead)~Gender+Smoke+Ctype+Cstage+Age_blood+radio+chemo+",metab)),
                       data=data_sub_st)
        cox.fit <- summary(model)
        lp <- predict(model,type="lp",newdata = data_sub_st)
        
        ROC <- timeROC(T=data_sub_st$OS_diag,
                       delta=data_sub_st$Dead,marker=lp,
                       cause=1,weighting="marginal",iid=T,
                       times=c(3,5,10))
      }else if(coltofilter[k]=="Gender"){
        model <- coxph(as.formula(paste0("Surv(OS_diag,Dead)~Smoke+Ctype+Cstage+Age_blood+radio+chemo+",metab)),
                       data=data_sub_st)
        cox.fit <- summary(model)
        lp <- predict(model,type="lp",newdata = data_sub_st)
        
        ROC <- timeROC(T=data_sub_st$OS_diag,
                       delta=data_sub_st$Dead,marker=lp,
                       cause=1,weighting="marginal",iid=T,
                       times=c(3,5,10))
      }else if(coltofilter[k]=="Smoke"){
        model <- coxph(as.formula(paste0("Surv(OS_diag,Dead)~Gender+Ctype+Cstage+Age_blood+radio+chemo+",metab)),
                       data=data_sub_st)
        cox.fit <- summary(model)
        lp <- predict(model,type="lp",newdata = data_sub_st)
        
        ROC <- timeROC(T=data_sub_st$OS_diag,
                       delta=data_sub_st$Dead,marker=lp,
                       cause=1,weighting="marginal",iid=T,
                       times=c(3,5,10))
      }else if(coltofilter[k]=="Ctype"){
        model <- coxph(as.formula(paste0("Surv(OS_diag,Dead)~Gender+Smoke+Cstage+Age_blood+radio+chemo+",metab)),
                       data=data_sub_st)
        cox.fit <- summary(model)
        lp <- predict(model,type="lp",newdata = data_sub_st)
        
        ROC <- timeROC(T=data_sub_st$OS_diag,
                       delta=data_sub_st$Dead,marker=lp,
                       cause=1,weighting="marginal",iid=T,
                       times=c(3,5,10))
      }else if(coltofilter[k]=="radio"){
        model <- coxph(as.formula(paste0("Surv(OS_diag,Dead)~Gender+Smoke+Cstage+Ctype+Age_blood+chemo+",metab)),
                       data=data_sub_st)
        cox.fit <- summary(model)
        lp <- predict(model,type="lp",newdata = data_sub_st)
        
        ROC <- timeROC(T=data_sub_st$OS_diag,
                       delta=data_sub_st$Dead,marker=lp,
                       cause=1,weighting="marginal",iid=T,
                       times=c(3,5,10))
      }else if(coltofilter[k]=="chemo"){
        model <- coxph(as.formula(paste0("Surv(OS_diag,Dead)~Gender+Smoke+Ctype+Cstage+Age_blood+radio+",metab)),
                       data=data_sub_st)
        cox.fit <- summary(model)
        lp <- predict(model,type="lp",newdata = data_sub_st)
        
        ROC <- timeROC(T=data_sub_st$OS_diag,
                       delta=data_sub_st$Dead,marker=lp,
                       cause=1,weighting="marginal",iid=T,
                       times=c(3,5,10))
      }else if(coltofilter[k]=="Class_stage"){
        model <- coxph(as.formula(paste0("Surv(OS_diag,Dead)~Gender+Smoke+Ctype+Cstage+Age_blood+radio+chemo+",metab)),
                       data=data_sub_st)
        cox.fit <- summary(model)
        lp <- predict(model,type="lp",newdata = data_sub_st)
        
        ROC <- timeROC(T=data_sub_st$OS_diag,
                       delta=data_sub_st$Dead,marker=lp,
                       cause=1,weighting="marginal",iid=T,
                       times=c(3,5,10))
      }
      
      #print(cox.fit)
      multi_betas[m,"beta"] <- cox.fit$coefficients[metab,1]
      multi_betas[m,"se"] <- cox.fit$coefficients[metab,3]
      multi_betas[m,"pvalue"] <- cox.fit$coefficients[metab,5]
      multi_betas[m,"z_stat"] <- cox.fit$coefficients[metab,4]
      multi_betas[m,"number"] <- nrow(data_sub_st)
      multi_betas[m,'AUC3'] <- ROC$AUC[1]
      multi_betas[m,'AUC5'] <- ROC$AUC[2]
      multi_betas[m,'AUC10'] <- ROC$AUC[3]
      
    }
    
    HRstrat <- multi_betas
    HRstrat <- HRstrat%>%mutate(Variable=coltofilter[k],strat_column=strat[i]) #add new columns for each strata 
    HR_strat_full <- bind_rows(HR_strat_full,HRstrat)
  }
}

HR_strat_full_main <- HR_strat_full[-1,]

HR_strat_full_main <- HR_strat_full_main%>%mutate(LCI=beta-1.96*se,UCI=beta+1.96*se)


library(dplyr)

strata <- c("Class_age","Gender","Smoke","Ctype","Class_stage",'radio','chemo')

strat_hetero_main <- chisq_strata_metab(HR_strat_full_main, metablist, strata)

tabs5 <- strat_hetero_main %>%
  select(Index, beta, LCI, UCI, pvalue, Variable, strat_column, number, phet) %>%
  rename(coltofilter = Variable, strat = strat_column) %>%
  left_join(subset(data_origin, select=c('Index','Compounds','Origin'))) %>%
  arrange(coltofilter,Origin,Compounds,strat) %>%
  mutate(coltofilter = recode(coltofilter,
                              "Class_age" = "Age",
                              "Class_stage" = "Stage",
                              "Ctype" = "Histology"),
         Compounds = case_when(
           Index == "MW0138609" ~ "K-3-(6\"Ac)Gal-7-Rha", 
           Index == 'MW0012126' ~ '11βOHAn',
           TRUE                 ~ Compounds
         )) %>%
  select(coltofilter,strat,number,Origin,Compounds,beta,LCI,UCI,pvalue,phet)

write.csv(tabs5, "03result/02Table/TabS5.csv", row.names = FALSE)

strat_hetero_main %>% filter(phet<=0.05) %>% View

#############Age########################################
data_sub <- metab_compute
data_sub$Class_age <- ifelse(data_sub$Age_blood<=65,'Young',ifelse(is.na(data_sub$Age_blood),NA,'Older'))
data_sub$Gender <- ifelse(data_sub$Gender==0,'Female',ifelse(is.na(data_sub$Gender),NA,'Male'))
data_sub$Smoke <- ifelse(data_sub$Smoke==0,'Never',
                         ifelse(is.na(data_sub$Smoke),NA,
                                ifelse(data_sub$Smoke==2,'Former','Current')))
data_sub$Class_stage <- ifelse(data_sub$Cstage %in% c(1,2),'Early',ifelse(data_sub$Cstage %in% c(3,4),'Late',NA))
data_sub$radio <- ifelse(data_sub$radio==0,'No',ifelse(is.na(data_sub$radio),NA,'Yes'))
data_sub$chemo <- ifelse(data_sub$chemo==0,'No',ifelse(is.na(data_sub$chemo),NA,'Yes'))

coltofilter <- c("Gender","Smoke","Ctype","Class_stage",'radio','chemo')

############ Age
metablist <- c(robust_exp_age$Index,robust_metab_age$Index)
HR_strat_full <- setNames(as.data.frame(matrix(nrow=1, ncol = 7)), c("Index","beta","se","pvalue","z_stat","Variable","strat_column"))
for (k in 1:length(coltofilter)) { 
  strat <- data_sub%>%pull(coltofilter[k])%>%as.factor()%>%levels 
  for (i in 1:length(strat)) { 
    data_sub_st <- data_sub%>%filter(!!as.symbol(coltofilter[k]) == strat[i]) 
    multi_betas <- data.frame(Index=NA,beta=NA, se=NA, pvalue=NA, z_stat=NA, number=NA, AUC3=NA,AUC5=NA,AUC10=NA)
    
    for (m in 1:length(metablist)){
      multi_betas[m,"Index"] <- metablist[m]
      metab <- metablist[m]
      #print(i)
      if (coltofilter[k]=="Class_age") {
        model <- coxph(as.formula(paste0("Surv(OS_diag,Dead)~Gender+Smoke+Ctype+Cstage+Age_blood+chemo+radio+Age_blood*",metab)),
                       data=data_sub_st)
        cox.fit <- summary(model)
        lp <- predict(model,type="lp",newdata = data_sub_st)
        
        ROC <- timeROC(T=data_sub_st$OS_diag,
                       delta=data_sub_st$Dead,marker=lp,
                       cause=1,weighting="marginal",iid=T,
                       times=c(3,5,10))
      }else if(coltofilter[k]=="Gender"){
        model <- coxph(as.formula(paste0("Surv(OS_diag,Dead)~Smoke+Ctype+Cstage+Age_blood+chemo+radio+Age_blood*",metab)),
                       data=data_sub_st)
        cox.fit <- summary(model)
        lp <- predict(model,type="lp",newdata = data_sub_st)
        
        ROC <- timeROC(T=data_sub_st$OS_diag,
                       delta=data_sub_st$Dead,marker=lp,
                       cause=1,weighting="marginal",iid=T,
                       times=c(3,5,10))
      }else if(coltofilter[k]=="Smoke"){
        model <- coxph(as.formula(paste0("Surv(OS_diag,Dead)~Gender+Ctype+Cstage+Age_blood+chemo+radio+Age_blood*",metab)),
                       data=data_sub_st)
        cox.fit <- summary(model)
        lp <- predict(model,type="lp",newdata = data_sub_st)
        
        ROC <- timeROC(T=data_sub_st$OS_diag,
                       delta=data_sub_st$Dead,marker=lp,
                       cause=1,weighting="marginal",iid=T,
                       times=c(3,5,10))
      }else if(coltofilter[k]=="Ctype"){
        model <- coxph(as.formula(paste0("Surv(OS_diag,Dead)~Gender+Smoke+Cstage+Age_blood+chemo+radio+Age_blood*",metab)),
                       data=data_sub_st)
        cox.fit <- summary(model)
        lp <- predict(model,type="lp",newdata = data_sub_st)
        
        ROC <- timeROC(T=data_sub_st$OS_diag,
                       delta=data_sub_st$Dead,marker=lp,
                       cause=1,weighting="marginal",iid=T,
                       times=c(3,5,10))
      }else if(coltofilter[k]=="radio"){
        model <- coxph(as.formula(paste0("Surv(OS_diag,Dead)~Gender+Smoke+Ctype+Cstage+Age_blood+chemo+Age_blood*",metab)),
                       data=data_sub_st)
        cox.fit <- summary(model)
        lp <- predict(model,type="lp",newdata = data_sub_st)
        
        ROC <- timeROC(T=data_sub_st$OS_diag,
                       delta=data_sub_st$Dead,marker=lp,
                       cause=1,weighting="marginal",iid=T,
                       times=c(3,5,10))
      }else if(coltofilter[k]=="chemo"){
        model <- coxph(as.formula(paste0("Surv(OS_diag,Dead)~Gender+Smoke+Ctype+Cstage+Age_blood+radio+Age_blood*",metab)),
                       data=data_sub_st)
        cox.fit <- summary(model)
        lp <- predict(model,type="lp",newdata = data_sub_st)
        
        ROC <- timeROC(T=data_sub_st$OS_diag,
                       delta=data_sub_st$Dead,marker=lp,
                       cause=1,weighting="marginal",iid=T,
                       times=c(3,5,10))
      }else if(coltofilter[k]=="Class_stage"){
        model <- coxph(as.formula(paste0("Surv(OS_diag,Dead)~Gender+Smoke+Ctype+Cstage+Age_blood+radio+chemo+Age_blood*",metab)),
                       data=data_sub_st)
        cox.fit <- summary(model)
        lp <- predict(model,type="lp",newdata = data_sub_st)
        
        ROC <- timeROC(T=data_sub_st$OS_diag,
                       delta=data_sub_st$Dead,marker=lp,
                       cause=1,weighting="marginal",iid=T,
                       times=c(3,5,10))
      }
      

      multi_betas[m,"beta"] <- cox.fit$coefficients[paste0('Age_blood:',metab),1]
      multi_betas[m,"se"] <- cox.fit$coefficients[paste0('Age_blood:',metab),3]
      multi_betas[m,"pvalue"] <- cox.fit$coefficients[paste0('Age_blood:',metab),5]
      multi_betas[m,"z_stat"] <- cox.fit$coefficients[paste0('Age_blood:',metab),4]
      multi_betas[m,"number"] <- nrow(data_sub_st)
      multi_betas[m,'AUC3'] <- ROC$AUC[1]
      multi_betas[m,'AUC5'] <- ROC$AUC[2]
      multi_betas[m,'AUC10'] <- ROC$AUC[3]
    }
    
    HRstrat <- multi_betas
    HRstrat <- HRstrat%>%mutate(Variable=coltofilter[k],strat_column=strat[i]) #add new columns for each strata 
    HR_strat_full <- bind_rows(HR_strat_full,HRstrat)
  }
}


HR_strat_full_age <- HR_strat_full[-1,]

HR_strat_full_age <- HR_strat_full_age%>%mutate(LCI=beta-1.96*se,UCI=beta+1.96*se)


library(dplyr)

strata <- c("Gender","Smoke","Ctype","Class_stage",'radio','chemo')

strat_hetero_age <- chisq_strata_metab(HR_strat_full_age, metablist, strata)

tabs6 <- strat_hetero_age %>%
  select(Index, beta, LCI, UCI, pvalue, Variable, strat_column, number, phet) %>%
  rename(coltofilter = Variable, strat = strat_column) %>%
  left_join(subset(data_origin, select=c('Index','Compounds','Origin'))) %>%
  arrange(coltofilter,Origin,Compounds,strat) %>%
  mutate(coltofilter = recode(coltofilter,
                              "Class_age" = "Age",
                              "Class_stage" = "Stage",
                              "Ctype" = "Histology"),
         ) %>%
  select(coltofilter,strat,number,Origin,Compounds,beta,LCI,UCI,pvalue,phet)

write.csv(tabs6, "03result/02Table/TabS6.csv", row.names = FALSE)                                         

strat_hetero_age %>% filter(phet<=0.05) %>% View #%>% distinct(protein, .keep_all = TRUE)

