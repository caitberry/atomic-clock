############################################################################
############################################################################
### Bayes
############################################################################
############################################################################

# stanDat = stan_data
# trueTau = tau
# trueMu = mu
# model = myHSSG_noDoF

#############################################################################################
#### Gauss + Gauss model
#############################################################################################

fitGaussGauss=function(stanDat,trueTau,trueMu,model,n.iter,adapt.delta,tree.depth){
  GGfit=sampling(model,
                 data=stanDat,
                 iter=n.iter, chains=3, cores=3,
                 verbose = F,refresh=0,control = list(adapt_delta = adapt.delta, max_treedepth = tree.depth))
  
  GGfit_df=as.data.frame(GGfit)

  mumcmc=GGfit_df$`mu`
  taumcmc=GGfit_df$`tau`
  
  tauout=data.frame(method="GaussGauss",
                    parameter="tau",
                    estimate=mean(taumcmc),
                    ci.lb=as.numeric(quantile(taumcmc,probs = .025)),
                    ci.ub=as.numeric(quantile(taumcmc,probs = .975)),
                    bad.rhat=sum(summary(GGfit)$summary[,"Rhat"]>1.05),
                    inInt=(as.numeric(quantile(taumcmc,probs = .025))<=trueTau &&
                             as.numeric(quantile(taumcmc,probs = .975)) >=trueTau),
                    min.neff=min(summary(GGfit)$summary[,"n_eff"]),
                    num.Divergent=get_num_divergent(GGfit),
                    tree.Depth.Ex=get_num_max_treedepth(GGfit),
                    min.bfmi=min(get_bfmi(GGfit)))
  muout=data.frame(method="GaussGauss",
                   parameter="mu",
                   estimate=mean(mumcmc),
                   ci.lb=as.numeric(quantile(mumcmc,probs = .025)),
                   ci.ub=as.numeric(quantile(mumcmc,probs = .975)),
                   inInt=(as.numeric(quantile(mumcmc,probs = .025))<=trueMu &&
                            as.numeric(quantile(mumcmc,probs = .975)) >=trueMu),
                   bad.rhat=sum(summary(GGfit)$summary[,"Rhat"]>1.05),
                   min.neff=min(summary(GGfit)$summary[,"n_eff"]),
                   num.Divergent=get_num_divergent(GGfit),
                   tree.Depth.Ex=get_num_max_treedepth(GGfit),
                   min.bfmi=min(get_bfmi(GGfit)))
  
  return(rbind(tauout,muout))
}
