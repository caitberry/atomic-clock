linearOP = function(x,u,nu=NA,m=1e4){
  #nu = numbers of degrees of freedom on which the standard uncertainties are based 
  weights = rep(1, length(x))

  index = sample(1:length(x), size=m, prob=weights/sum(weights), replace=TRUE)
  mu = x[index]
  sigma = u[index]
  if (!anyNA(nu) && length(nu)!=0){

    nuSample = nu[index]

    factorVar = numeric(m)
    y=numeric(m)

    indexDFgr2 = (nuSample > 2 & nuSample !=Inf)
    indexDFless2 = (nuSample <= 2)

    indexInf = (nuSample == Inf)


    phi=numeric(m)
    a = 0.4668994
    b = -0.3998882
    phi = 1.5*((1 - (3/4)*(a - 4*b *(nuSample^(-3/4)-1)/3))^(-4/3))


    factorVar[indexDFgr2] = sigma[indexDFgr2]/
      sqrt(nuSample[indexDFgr2]/(nuSample[indexDFgr2]-2))

    factorVar[indexDFless2] = sigma[indexDFless2]/phi[indexDFless2]

    y[indexDFgr2] = mu[indexDFgr2] + rt(sum(indexDFgr2), df=nuSample[indexDFgr2]) * factorVar[indexDFgr2]

    y[indexDFless2] = mu[indexDFless2] + rt(sum(indexDFless2), df=nuSample[indexDFless2]) *  factorVar[indexDFless2]

    y[indexInf]=rnorm(sum(indexInf), mean=mu[indexInf], sd=sigma[indexInf])

  }else{
    y = rnorm(m, mean=mu, sd=sigma)

  }

  return(y)
}

#mean(res) #consensus value
#sd(res) #standard uncertainty
#quantile(res, c(0.025, 0.975)) #95% coverage interval