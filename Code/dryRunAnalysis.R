######################################################################
### Read in data 
######################################################################
rm(list=ls())

library(readr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(openxlsx)
library(fields)

library(tidyverse) #data wrangling
library(RSpectra) #eigenvalue solving

mtse=modules::use("Functions.R")

###############################################################################################
###############################################################################################
dat2024_04_24_ErYb_AlSr <- read.csv("Data/2024-04-24-ErYb-AlSr.csv",header = F,sep="",na.strings = "NaN")
dim(dat2024_04_24_ErYb_AlSr)
colnames(dat2024_04_24_ErYb_AlSr)=c("MJD","FracDiff")
dat2024_04_24_ErYb_AlSr$seconds=1:dim(dat2024_04_24_ErYb_AlSr)[1]


# source("Functions.R")
dat=dat2024_04_24_ErYb_AlSr%>%mutate(missing=is.na(FracDiff),
                                     date=convertToDateTime(MJD, origin = "1858-11-17",tz="MDT"))


dat=filter(dat, missing==F)
t.vec <- dat$seconds
x.t=dat$FracDiff
N=length(x.t)

test=spectralEstWithUnc(x.t = x.t,t.vec=t.vec,N.fourier = floor(N/2) + 1,#100,
                        numTapers = 10,calcCov = F,
                        myW = 4/N*3)
test$V.mat$e.values
V.mat <- get_tapers(t.vec, W = 4/N*3, K = 9) 
V.mat$e.values



### data analysis steps
# 1. read in data
# 2. concatenate small gaps
# 3. get spectral estimate
# 3a) try a variety of W and K values and check eigenvalues to make sure you've made good selections
# can start with 8/N or 12/N for W
# 4. Look at data, spectral estimate, and tapers
# 5. calculate avar estimate with spectrum for a series of tau values

# 6. calculate avar estimate using old method for same tau series
# 7. plot both with uncertainties

# 
# 
# 
# 
# #############################################old below
# 
# 
# 
# #####################make a function
# source("Functions.R")
# 
# calcAvar=function(thetau,freq,spec.hat,delta.f){
#   # G_tau vector length number of frequencies 
#   G_tau <- transfer.func(freq,tau = thetau) #change the tau value to get different vectors
#   G_tau[1] <- 0 # this was 1 in the old code, but should be 0
#   
#   avar=G_tau%*%spec.hat*delta.f
#   
#   return(avar)
# }
# 
# calculateAvars=function(x.t,t.vec,taus,N.fourier=100,numTapers=3,calcCov,myW){
#   N <- length(t.vec)
#   # N.fourier <- 1000#floor(N/2) + 1
#   freq <- seq(0,0.5, length.out = N.fourier)
#   
#   delta.f <- freq[2] #interval spacing between frequencies, needed for spectral avar calculation
#   
#   ##calculate tapers for this data spacing
#   V.mat <- get_tapers(t.vec, W = myW, K = numTapers) #W was 4/N
#   
#   MTSE_full <- MT_spectralEstimate_freqs(x.t, freq, V.mat$tapers) 
#   
#   ### calculate the covariance matrix 
#   Cov.mat_chave <- matrix(NA, nrow = N.fourier, ncol = N.fourier)
#   
#   if(calcCov==T){
#     for(i in 1:N.fourier){
#       if(i %% 100 == 0){print(paste(i," of ",N.fourier))}
#       j = 1
#       while(j <= i){
#         Cov.mat_chave[i,j] <- norm(Conj(t(V.mat$tapers*exp(-im*2*pi*freq[i]*t.vec)*(1/sqrt(numTapers))))%*%(V.mat$tapers*exp(-im*2*pi*freq[j]*t.vec)*(1/sqrt(numTapers))), type = "2") 
#         j = j+1
#       }
#     }
#     
#     Cov.mat_chave[upper.tri(Cov.mat_chave)] <- t(Cov.mat_chave)[upper.tri(Cov.mat_chave)]
#   }
#   
#   
#   ### calc avar and variance of avar
#   
#   cov.mat=avar=numeric(length(taus))
#   
#   for(i in 1:length(taus)){
#     #taken from above, need to calculate the covariance
#     G_tau <- transfer.func(freq,tau = taus[i]) 
#     G_tau[1] <- 0 
#     
#     avar[i]=calcAvar(taus[i],MTSE_full$freqs,MTSE_full$spectrum,delta.f)
#     
#     #calculate variance for the AVAR estimate at the given tau
#     cov.mat[i] <- t(G_tau)%*%(Cov.mat_chave)%*%G_tau*(delta.f)^2#/(sqrt(numTapers)) #can't tell if this term should be here, details below
#   }
#   
#   avarOut=data.frame(tau=taus,avar=avar,var=cov.mat)
#   
#   ############### analyze old way
#   avar.calc <- getAvars(N,x.t, taus = taus)
#   oldRes=data.frame(tau=avar.calc$avarRes$taus,avar=avar.calc$avarRes$avars,calculation="old")
#   overRes=data.frame(tau=avar.calc$avarRes$taus,avar=avar.calc$avarRes$overavars,calculation="overlapping")
#   spectralRes=data.frame(tau=avar.calc$avarRes$taus,avar=avarOut$avar,calculation="spectral")
#   allRes=bind_rows(oldRes,overRes,spectralRes)
#   
#   return(list(x.t=x.t,t.vec=t.vec,V.mat=V.mat,MTSE_full=MTSE_full,
#               freq=freq,Cov.mat_chave=Cov.mat_chave,avarOut=avarOut,allAvarRes=allRes,W=myW))
# }
