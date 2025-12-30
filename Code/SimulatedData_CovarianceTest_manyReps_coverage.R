library(tidyr)
library(ggplot2)
library(dplyr)

numberOfSimulations = 100

set.seed(123)
# N=2048+2048*.5 #length of data with gaps
N=200

# omitted=c(100:400, 750:864 ,1300:1600, 1700:1905, 2013:2113)#1024 missing
omitted=c(100:140)

t.vec <- 1:N
t.vec[omitted] <- NA

mtse=modules::use("Functions.R")
source("Functions_SpectrumCovariance.R")

## saving the date to label file outputs
runDate=format(Sys.Date(),"%m%d%y")

setWnum = 12
setW = setWnum/N
setK = 5
V.mat <- mtse$get_tapers(t.vec, W = setW, K = setK) 
print(V.mat$e.values)
taperMatrix=V.mat$tapers

taus <- 2^(0:9)
taus <- taus[taus<floor(N/3)]

avarDF=data.frame()


N.omitted=length(na.omit(t.vec))
# edfOLD=c()
# i=1
# for(m in taus){
#   edfOLD[i] <- ((3*(N.omitted-1)/(2*m)) - (2*(N.omitted-2)/N.omitted))*(4*m^2)/(4*m^2 + 5)
#   i=i+1
# }

for(i in 1:numberOfSimulations){
  print(i)
  x.t <- rnorm(N)  
  x.t[omitted]=NA
  
  MTSE_full <- mtse$MT_spectralEstimate_fft(x.t, taperMatrix) 
  specAVARest=mtse$AVAR_spec(spectral_est = MTSE_full$spectrum,taus = taus,calcUnc = F)
  
  # Cov.mat=spec_cov.mat(x.t = x.t, t.vec = t.vec, N.fourier = floor(N/2) + 1, 
  #                      taperMat = V.mat$tapers, isWhite = F,
  #                      max.lag.acf = 4,
  #                      numCores=4)
  
  Cov.matWN=spec_cov.mat_slow_WN(t.vec = t.vec, N.fourier = floor(N/2) + 1, 
                       taperMat = V.mat$tapers)
  
  # specAVARest=mtse$AVAR_spec(spectral_est = MTSE_full$spectrum,taus = taus,calcUnc = T,Cov.mat = Cov.mat)
  specAVARest=mtse$AVAR_spec(spectral_est = MTSE_full$spectrum,taus = taus,calcUnc = T,Cov.mat = Cov.matWN)
  
  oldAVARest=mtse$overlapping_avar_fn(y = na.omit(x.t),m = taus)
  oldAVARestUncertainty=mtse$avar_CI(CI.level = .68,
                                     noise_type = "white noise",
                                     avar_type = "chisquared",
                                     avars = oldAVARest,
                                     taus=taus,
                                     N=N.omitted)
  # #########
  # spectral_edf=edfOLD
  spectral_edf=2*(specAVARest$avar)^2/specAVARest$avarVar
  # a=0.6826895
  # qnorm((1-a)/2)
  # qnorm((1+a)/2)
  CI.level=.68
  a <- (1-CI.level)/2

  spectral_CI.limits <- dplyr::bind_rows("lower" = specAVARest$avar*spectral_edf/stats::qchisq(1-a,spectral_edf),
                                         "upper" = specAVARest$avar*spectral_edf/stats::qchisq(a, spectral_edf))
  

  avarDFnew=data.frame(avar=c(oldAVARest,specAVARest$avar),
                    avarLower=c(oldAVARestUncertainty$lower,spectral_CI.limits$lower),
                    avarUpper=c(oldAVARestUncertainty$upper,spectral_CI.limits$upper),
                    tau=rep(taus,2),
                    Method=rep(c("Current","Spectral"),each=length(taus)),
                    rep=i)
  avarDF=bind_rows(avarDF,avarDFnew)

}

saveRDS(avarDF,file = "Results/avarDFforcoverageProb012125WN.Rds")



################ read and plot


lotsOfAvarReps=readRDS("Results/avarDFforcoverageProb011725.Rds")
head(lotsOfAvarReps)

lotsOfAvarReps=lotsOfAvarReps %>% group_by(tau, Method) %>%
  mutate(truth=1/tau,coverage=(avarLower < truth & avarUpper > truth))
View(lotsOfAvarReps)

numberOfSimulations=max(lotsOfAvarReps$rep)

ggplot(lotsOfAvarReps %>% summarise(covProb=sum(coverage)/numberOfSimulations),
       aes(tau,covProb,col=Method))+
  geom_point()+
  geom_hline(yintercept = .68)


ggplot(filter(lotsOfAvarReps,tau==32,rep<100),aes(rep,avar,ymin=avarLower,ymax=avarUpper,color=coverage))+
  geom_point()+
  geom_errorbar()+
  facet_wrap(~Method)


data=1/filter(lotsOfAvarReps,tau==1,Method=="Spectral")$avar
# data=rchisq(1000,df = 10,ncp = 0)

var(data)
hist(data)
# estimate_chisq_params <- function(data) {
  # Negative log-likelihood function
  negLogLik <- function(df, ncp) {
    -sum(dchisq(data, df = df, ncp = ncp, log = TRUE))
  }
  
  # Initial values (adjust based on your data if necessary)
  start_vals <- list(df = trunc(var(data)/2), ncp = trunc(mean(data)))
  # start_vals <- list(df = 10, ncp = 1)
  
  # MLE using bbmle
  mle_fit <- bbmle::mle2(negLogLik, start = start_vals)
  mle_fit <- optim(start_vals,negLogLik)
  # Return estimated parameters as a named vector
  # df <- dplyr::tibble(
  #   est_df = coef(mle_fit)[1],
  #   est_ncp = coef(mle_fit)[2]
  # )
  # return(df)
# }

estimate_chisq_params(test)
