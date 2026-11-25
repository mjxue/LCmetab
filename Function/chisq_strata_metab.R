
chisq_strata_metab <- function(data,metablist, strata){
  data2 <- data%>%mutate(phet=NA)%>%mutate(phet=as.numeric(phet))
  for (i in 1:length(metablist)) {
    df <- data%>%filter(Index==metablist[i])
    for(j in 1:length(strata)){
      data_strat <- df%>%filter(Variable==strata[j])%>%filter(!is.na(beta)&!is.na(se))
      yi <- data_strat$beta # observed effect size estimates
      vi <- data_strat$se 
      vi <- vi^2 # corresponding sampling variances
      het <- tryCatch(metafor::rma(yi, vi, method="FE"), error=function(e){return(NA)}) # fixed model effect
      if(!is.na(het$QEp)){
        p_het <- het$QEp
      }else{
        p_het <- NA
      }
      data2[data2$Index==metablist[i]&data2$Variable==strata[j],"phet"] <- p_het
    }
    
  }
  return(data2)
}