######################################################################
##
## FILE        : modifiedKHmethod.R

## AUTHOR      : Amanda Koepke & Antonio Possolo
## MODIFICATION: 2019-Apr-12

## INPUT : x   = Numeric vector with values measured by the labs (required)
##         u   = Numeric vector with standard uncertainties
##               associated with measured values (required)
##         coverageProb = Coverage probability for the interval
##               (DoE.Lwr, DoE.Upr), and for the interval DoE.x +/-
##               DoE.U (required)

## OUTPUT: Results of DL with modified
##         version of the Hartung-Knapp-Sidik-Jonkman
##         approximation, as described by Rover, Knapp, 
##         and Friede (2015, DOI 10.1186/s12874-015-0091-1)

######################################################################

modifiedKHmethod = function (x, u, coverageProb) 
{
  require(metafor)
  
  n = length(x)
  w0 = 1/u^2
  x0 = sum(w0*x)/sum(w0)
  tau2 = max(0, (sum(w0*(x-x0)^2)-n+1) /
               (sum(w0) - sum(w0^2)/sum(w0)))
  w = 1/(u^2+tau2)
  mu = sum(w*x)/sum(w)
  # q=sum(w*(x-mu)^2)/(n-1) # adjustment calculated the old way
  qSTAR = max(1, sum(w*(x-mu)^2)/(n-1)) # max(1,q)
  mu.u = sqrt(qSTAR/sum(w))
  mu.df = n-1
  
  ci.ub=mu+(sqrt(qSTAR)*mu.u*qt((1+coverageProb)/2,mu.df)) ## ci.ub for modified adjustment
  ci.lb=mu+(sqrt(qSTAR)*mu.u*qt((1-coverageProb)/2,mu.df)) ## ci.lb for modified adjustment
  
  DLout=data.frame(b=mu,
                   se=mu.u,
                   ci.ub=ci.ub,
                   ci.lb=ci.lb,
                   tau2=tau2)
  return(DLout)
}