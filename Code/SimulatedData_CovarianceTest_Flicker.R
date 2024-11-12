library(RobPer)

set.seed(123)
N=1000 #length of data with gaps

# x.t <- rnorm(N) + 0.001*TK95(N = N, alpha = 1) 
x.t <- TK95(N = N, alpha = 1)
omitted=c(301:400)
x.t[omitted]=NA

t.vec <- 1:N
t.vec[which(is.na(x.t))] <- NA

mtse=modules::use("Functions.R")
source("Functions_SpectrumCovariance.R")

W=4/N*3
K=5
V.mat <- mtse$get_tapers(t.vec, W = W, K = K) 
print(V.mat$e.values)
taperMatrix=V.mat$tapers

MTSE_full <- mtse$MT_spectralEstimate_fft(x.t, taperMatrix) 
plot(t.vec,x.t)
plot(MTSE_full$freqs,MTSE_full$spectrum)

Cov.mat=spec_cov.mat(x.t = x.t, t.vec = t.vec, N.fourier = floor(N/2) + 1, 
                     taperMat = V.mat$tapers, isWhite = F,
                     max.lag.acf = 4,
                     numCores=10)
diag(Cov.mat)[1]
specDF=data.frame(frequency=MTSE_full$freqs,
                  S=MTSE_full$spectrum,
                  unc=sqrt(diag(Cov.mat)[1]))
library(ggplot2)
ggplot(specDF,aes(frequency,S,ymin=S-2*unc,ymax=S+2*unc))+
  geom_point()+
  geom_errorbar()

# 5. calculate avar estimate with spectrum for a series of tau values


taus <- 2^(0:9)
taus <- taus[taus<floor(N/3)]

specAVARest=mtse$AVAR_spec(spectral_est = MTSE_full$spectrum,taus = taus,calcUnc = T,Cov.mat = Cov.mat)

# 6. calculate avar estimate using old method for same tau series

plot(x.t)

oldAVARest=mtse$overlapping_avar_fn(y = na.omit(x.t),m = taus)

oldAVARestUncertainty=mtse$avar_CI(CI.level = .68,
                                   noise_type = "white noise", 
                                   avar_type = "chisquared", 
                                   avars = oldAVARest, 
                                   taus=taus,
                                   N=N)

# 7. plot both with uncertainties


avarDF=data.frame(avar=c(oldAVARest,specAVARest$avar),
                  avarLower=c(oldAVARestUncertainty$lower,specAVARest$avar-sqrt(specAVARest$avarVar)),
                  avarUpper=c(oldAVARestUncertainty$upper,specAVARest$avar+sqrt(specAVARest$avarVar)),
                  tau=rep(taus,2),
                  Method=rep(c("Current","Spectral"),each=length(taus)))

ggplot(avarDF,aes(tau,avar,col=Method,ymin=avarLower,ymax=avarUpper))+
  geom_point()+
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  # theme(legend.position = c(.15, .2))+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()+
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))+
  geom_errorbar()
# facet_wrap(~Ratio)

library(dplyr)
ggplot(filter(avarDF,Method=='Spectral'),aes(tau,avar,col=Method,ymin=avarLower,ymax=avarUpper))+
  geom_point()+
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  # theme(legend.position = c(.15, .2))+
  scale_y_log10()+
  scale_x_log10()+
  annotation_logticks()+
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))+
  geom_errorbar()
# facet_wrap(~Ratio)

