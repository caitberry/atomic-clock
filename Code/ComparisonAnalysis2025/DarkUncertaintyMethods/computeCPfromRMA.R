analyzeDataRMA=function(data,trueTau,trueMu){
  ## more on CI for tau produced by metafor
  # https://stats.stackexchange.com/questions/167497/the-confidence-interval-of-tau2-and-i2-produced-by-the-meta-and-metafor-p
  
  ## CI may not include estimate. Should I use a different way to get the CI? Seems to only be a problem for DL 
  # https://stats.stackexchange.com/questions/476266/unusual-problem-with-metaforrma-output
  
  DLout=rma(yi = data$x,sei = data$u,method = "DL")
  SJout=rma(yi = data$x,sei = data$u,method = "SJ")
  REMLout=rma(yi = data$x,sei = data$u,method = "REML")
  PMout=rma(yi = data$x,sei = data$u,method = "PM",control=list(tau2.max=1000))
  
  CIs=rbind(data.frame(as.list(confint(DLout)$random["tau",])),
            data.frame(as.list(confint(SJout)$random["tau",])),
            data.frame(as.list(confint(REMLout)$random["tau",])),
            data.frame(as.list(confint(PMout)$random["tau",])))
  
  CImethods=data.frame(method=c("DL","SJ","REML","PM"),
                       parameter=c("tau","tau","tau","tau"))
  
  allTauCIout=cbind(CImethods,CIs %>% mutate(inInt=ci.lb<=trueTau && ci.ub>=trueTau))
  
  mu.est=c(DLout$b, SJout$b, REMLout$b, PMout$b)
  mu.ci.lb=c(DLout$ci.lb, SJout$ci.lb, REMLout$ci.lb, PMout$ci.lb)
  mu.ci.ub=c(DLout$ci.ub, SJout$ci.ub, REMLout$ci.ub, PMout$ci.ub)
  
  muCI = data.frame(method=c("DL","SJ","REML","PM"),
                    parameter="mu",
                    estimate=mu.est,
                    ci.lb=mu.ci.lb,
                    ci.ub=mu.ci.ub)
  allMuCIout=cbind(muCI %>% mutate(inInt=ci.lb<=trueMu && ci.ub>=trueMu))
  
  bothCIout=rbind(allTauCIout,allMuCIout)
  placeholdersForBayesOut=data.frame(bad.rhat=0,
                                     min.neff=NA,
                                     num.Divergent=NA,
                                     tree.Depth.Ex=NA,
                                     min.bfmi=NA)
  return(cbind(bothCIout,placeholdersForBayesOut))
}
