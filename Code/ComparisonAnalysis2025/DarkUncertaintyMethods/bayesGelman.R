
bayesGelman = function(x, u, nu, 
                       tauPriorScale, sigmaPriorScale,
                       ni, nb, nt, nc=1){
  require(R2jags)
  n = length(x)
  
  
  if (is.null(nu)) {
    ## When nu is NULL we take sigma=u, as if all
    ## elements of u were based on infinitely many
    ## degrees of freedom
    cat("model {
        for(j in 1:n)
        {
        xi[j] ~ dnorm(mu, tauINV2)
        uINV2[j] <- pow(u[j],-2)
        x[j] ~ dnorm(xi[j], uINV2[j])
        }
        mu ~ dnorm(0.0, 1.0E-10)
        tau.prvar <- pow(tauPriorScale, -2)
        tau ~ dt(0, tau.prvar, 1)T(0,)
        tauINV2 <- pow(tau,-2)
  }", fill=TRUE, file="LabComp-Model.txt")
    LabComp.Data = list(n=n, x=x, u=u, tauPriorScale=tauPriorScale)
    
    ## Inits function
    LabComp.inits = function() {
      list(mu = rnorm(1, mean=median(x), sd=mad(x)),
           xi = rnorm(n, mean=median(x), sd=u) ) }
    
    ## Parameters to estimate
    LabComp.params = c("mu", "tau", "xi")
  } else {
    ## nu that are Inf, then we replace Inf by 1000
    nu[nu == Inf] = 1000
    cat("model {
      for(j in 1:n)
      {
      xi[j] ~ dnorm(mu, tauINV2)
      sigma[j] ~ dt(0, sig.prvar, 1)T(0,)
      sigmaINV2[j] <- pow(sigma[j],-2)
      x[j] ~ dnorm(xi[j], sigmaINV2[j])
      varmod[j]~ dgamma(nu[j]/2, 1/(2*sigma[j]*sigma[j]))
      }
      sig.prvar<-pow(sigmaPriorScale,-2)
      mu ~ dnorm(0.0, 1.0E-10)
      tau.prvar <- pow(tauPriorScale,-2)
      tau ~ dt(0, tau.prvar, 1)T(0,)
      tauINV2 <- pow(tau,-2)
}", fill=TRUE, file="LabComp-Model.txt")
    LabComp.Data = list(n=n, x=x, nu=nu, varmod=nu*u^2,
                        tauPriorScale=tauPriorScale, sigmaPriorScale=sigmaPriorScale )
    
    ## Inits function
    LabComp.inits = function() {
      list(mu = rnorm(1, mean=median(x), sd=mad(x)),
           xi = rnorm(n, mean=median(x), sd=u)) }
    
    ## Parameters to estimate
    LabComp.params = c("mu", "tau", "xi", "sigma")
  }
  
  ## Start Gibbs sampling
  
  withProgress(message = 'Running Bayes', style="old",value = 0, {
    LabComp.output = jags(LabComp.Data,
                          inits=LabComp.inits, parameters.to.save=LabComp.params,
                          model.file="LabComp-Model.txt", n.thin=nt, n.chains=nc,
                          n.burnin=nb, n.iter=ni)
  })
  bayesres = as.mcmc(LabComp.output)
  geweke = geweke.diag(bayesres)
  zScores = geweke[[1]][[1]][-1] # confirm that this index is correct 
  # for (jc in 2:nc) {zScores = c(zScores, geweke[[jc]][[1]][-1])} #changed nc to 1
  geweke.pvalues = p.adjust(2*pnorm(-abs(zScores)), method="BH")
  if (min(geweke.pvalues) < 0.05) {
    cat("WARNING: MCMC may not have reached equilibrium\n")
    cat("WARNING: Results are unreliable\n")
    cat(paste("SUGGESTION: Re-run with ",
              "'Total number of iterations' = ", 10*ni, ", 'Length of burn in' = ", 2*nb,
              ", and 'Thinning rate' = ", ceiling(0.001*ni), "\n", sep="")) 
    
    # testWarn=paste("WARNING: MCMC may not have reached equilibrium\n","WARNING: Results are unreliable\n",paste("SUGGESTION: Re-run bayesGelman with ",
    #           "ni = ", 10*ni, ", nb = ", 2*nb,
    #           ", and nt = ", ceiling(0.001*ni), "\n", sep=""),sep=" ") 
    testWarn=paste("WARNING: MCMC may not have reached equilibrium<br/>",
                   "WARNING: Results are unreliable<br/>",
                   paste("SUGGESTION: Re-run with ",
                         "'Total number of iterations' = ", 10*ni, ", 'Length of burn in' = ", 2*nb,
                         ", and 'Thinning rate' = ", ceiling(0.001*ni), "<br/>", sep=""),sep=" ") 
    
  }else{
    testWarn=""
  }
  
  if (file.exists("LabComp-Model.txt")) {
    file.remove("LabComp-Model.txt")}
  
  s = summary(bayesres)
  
  return(list(mcmcout=bayesres,warn=testWarn))
  
}


