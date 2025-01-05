library(tidyr)
library(ggplot2)

numberOfSimulations = 100

set.seed(123)
N=2048+2048*.5 #length of data with gaps

omitted=c(100:400, 750:864 ,1300:1600, 1700:1905, 2013:2113)#1024 missing

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

samat=oamat=matrix(NA, nrow = numberOfSimulations, ncol = length(taus))

for(i in 1:numberOfSimulations){
  print(i)
  x.t <- rnorm(N)  # Replace with your time series data
  x.t[omitted]=NA
  
  MTSE_full <- mtse$MT_spectralEstimate_fft(x.t, taperMatrix) 
  specAVARest=mtse$AVAR_spec(spectral_est = MTSE_full$spectrum,taus = taus,calcUnc = F)
  samat[i,]=specAVARest$avar
  
  oldAVARest=mtse$overlapping_avar_fn(y = na.omit(x.t),m = taus)
  oamat[i,]=oldAVARest
  
}


saveRDS(samat,paste("Results/samat",runDate,"_W",setWnum,"_K",setK,"_N",N,"_",numberOfSimulations,"sims_WhiteNoiseGaps.Rds",sep=""))
saveRDS(oamat,paste("Results/oamat",runDate,"_W",setWnum,"_K",setK,"_N",N,"_",numberOfSimulations,"sims_WhiteNoiseGaps.Rds",sep=""))





##############################################################################
### Plotting just the sim results
### Doesn't yet look like Cait's final plots in style
##############################################################################

##tidy the data
both.res.mat <- rbind(oamat,samat)
method_labels <- rep(c("avar", "spectral"), each = numberOfSimulations)
df.messy <- as.data.frame(cbind(method_labels, both.res.mat))
colnames(df.messy) <-c("method", taus)


dat <- df.messy %>% gather(tau, measurement, -method)

dat$measurement <- as.numeric(dat$measurement)
dat$tau <- as.numeric(dat$tau)

breaks <- 10^(-10:10)
minor_breaks <- rep(1:9,21 )*rep(10^(-10:10), each = 9)

linedat_w=data.frame(tau=1:1000)
linedat_w$truth=(1/linedat_w$tau)#+linedat$tau
linedat_w$method="Truth"

ggplot(dat,aes(tau,measurement,col=method,group=interaction(tau,method)))+
  geom_boxplot(lwd = 1.2)+
  ### add true straight line below
  # geom_abline(slope = -1,intercept = 0,size=1)+
  ### add true curved line below, calculate beforehand!
  geom_line(data=linedat_w,aes(tau,truth), linewidth = 1)+
  ### This cahnges the legend title and labels
  scale_color_discrete(labels= c("avar"= "Avar","spectral"= "Spectral"),name="Method")+
  # ### all this makes the plot look more like a base R plot
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ### where to put legend. one option is "bottom", avoid having to place it. The tuple I have here
  ### basically specifies the x and y position in terms of the plot size in unit scale.
  theme_bw(base_size = 20)+
  theme(legend.position = c(.15, .2))+
  scale_y_log10(breaks = breaks, minor_breaks = minor_breaks)+
  scale_x_log10(breaks = breaks, minor_breaks = minor_breaks)+
  annotation_logticks()+
  ylab(expression(sigma^2*(tau)))+
  xlab(expression(tau))







##############################################################################
### Add error bars for one run
##############################################################################

#will just use the last run
plot(MTSE_full$freqs,MTSE_full$spectrum)
Cov.mat=spec_cov.mat(x.t = x.t, t.vec = t.vec, N.fourier = floor(N/2) + 1, 
                     taperMat = V.mat$tapers, isWhite = F,
                     max.lag.acf = 4,
                     numCores=10)
diag(Cov.mat)[1]
specAVARest=mtse$AVAR_spec(spectral_est = MTSE_full$spectrum,taus = taus,calcUnc = T,Cov.mat = Cov.mat)
# 
# # 6. calculate avar estimate using old method for same tau series
# 
# plot(x.t)
# 
# oldAVARest=mtse$overlapping_avar_fn(y = na.omit(x.t),m = taus)
# 
# oldAVARestUncertainty=mtse$avar_CI(CI.level = .68,
#                                    noise_type = "white noise", 
#                                    avar_type = "chisquared", 
#                                    avars = oldAVARest, 
#                                    taus=taus,
#                                    N=N)
# 
# # 7. plot both with uncertainties
# 
# 
# avarDF=data.frame(avar=c(oldAVARest,specAVARest$avar),
#                   avarLower=c(oldAVARestUncertainty$lower,specAVARest$avar-sqrt(specAVARest$avarVar)),
#                   avarUpper=c(oldAVARestUncertainty$upper,specAVARest$avar+sqrt(specAVARest$avarVar)),
#                   tau=rep(taus,2),
#                   Method=rep(c("Current","Spectral"),each=length(taus)))
# 
# ggplot(avarDF,aes(tau,avar,col=Method,ymin=avarLower,ymax=avarUpper))+
#   geom_point()+
#   ### add true straight line below
#   # geom_abline(slope = -1,intercept = 0,size=1)+
#   # theme(legend.position = c(.15, .2))+
#   scale_y_log10()+
#   scale_x_log10()+
#   annotation_logticks()+
#   ylab(expression(sigma^2*(tau)))+
#   xlab(expression(tau))+
#   geom_errorbar()
# # facet_wrap(~Ratio)
# 
