# source(file = "/home/aak3/NIST/atomic-clock/Code/ComparisonAnalysis2025/DarkUncertaintyMethods/coverageProbExp_withParallel_normal.R")
# .libPaths("/home/aak3/R/x86_64-redhat-linux-gnu-library/3.6/")

rsession=1

library(dplyr)

mypath="/home/aak3/NIST/atomic-clock/Code/ComparisonAnalysis2025/DarkUncertaintyMethods/"

############################################################################
### set parameters based on data
############################################################################

# df_AlSr <- read_csv("Data/ClockComparison2025/BayesianAnalysisData/ErYb_AlSr_data.csv")
# df_AlYb <- read_csv("Data/ClockComparison2025/BayesianAnalysisData/ErYb_AlYb_data.csv")
df_YbSr <- read_csv("Data/ClockComparison2025/BayesianAnalysisData/ErYb_YbSr_data.csv")

measurements <- df_YbSr$offset#c(1.000005, 1.000012, 0.999990, 1.000025, 0.999980)
uncertainties <- df_YbSr$statistical_unc#c(0.000005, 0.000005, 0.000004, 0.000006, 0.000005)

N <- length(measurements)

######

numberOfSims_set=10 ### make bigger

N_set=c(N,N+20)
tau_set=c(3,10)

allparams=expand.grid(N_set=N_set,tau_set=tau_set) #want every combo of N and tau

allparams$lb=min(uncertainties)
allparams$ub=max(uncertainties)

myparams=allparams %>%
  mutate(adaptDeltaVals = #########right now, setting these based on N
           ifelse(N_set == 5, .99,
                  ifelse(N_set== 10, .95,
                         ifelse(N_set== 20, .8,
                                ifelse(N_set== 50, .8,
                                       NA))))) %>%
  mutate(maxTreeDepthVals = #########right now, setting these based on N
           ifelse(N_set== 5, 17, 
                  ifelse(N_set== 10, 15,
                         ifelse(N_set== 20, 10,
                                ifelse(N_set== 50, 10,
                                       NA))))) %>%
  mutate(tdf_set =
           ifelse(N_set == 5 & tau_set %in%c(1,2), 2,
                  4)) #set prior sd to 4 for all except small N and small tau case, then too big


##############################################################################################
### don't change this (often)
##############################################################################################

mu_set=mean(measurements)
num.iters_set=5000
alpha_set=NA
nu_set=NA
howDatSim_set="normal"

rundate=format(Sys.Date(),"%y%m%d")

##############################################################################################
##############################################################################################
### simulate data
##############################################################################################
##############################################################################################

source(file = paste(mypath,"simulateData.R",sep=""))

############################################################################
############################################################################

# set.seed(4)
library("metafor")
library("dplyr")
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

######################  Analyze data

############################################################################
### analyze data from metafor
############################################################################

source(file = paste(mypath,"computeCPfromRMA.R",sep = ""))

############################################################################
### analyze data from Bayesian
############################################################################

source(file = paste(mypath,"fitBayesianModels.R",sep = ""))

############################################################################
### analyze data from new method
############################################################################

source(file = paste(mypath,"computeForCSmethods.R",sep = ""))

############################################################################
############################################################################
### do for many
############################################################################
############################################################################
source(file = paste(mypath,"runAllModels.R",sep = "")) ## will need to move this earlier and check for errors on the first run now

# numberOfSims=numberOfSims_set
# N=N_set[1]
# tau=tau_set[1]
# mu=mu_set
# num.iters=num.iters_set
# alpha=alpha_set
# nu=nu_set
# howDatSim=howDatSim_set
# rsession=1
# rsessionNum=1
# theadaptdelta=.99
# thetreedepth=15
# thetdf=4
# today=rundate
# lb=.1
# ub=1
# 



runForManyRepeats=function(numberOfSims,N,tau,mu,num.iters,alpha,nu,howDatSim,
                           rsessionNum,theadaptdelta,thetreedepth,thetdf,today,lb,ub){

  startTime=Sys.time()

  print(paste("tau is ",tau, ", N is ",N,", adapt delta is ",theadaptdelta,
              ", max tree is ",thetreedepth,", ub is ",ub,", lb is ",lb))

  myHGG_noDoF<-stan_model("/home/aak3/NIST/BTF/BTF_FY21/modelCode/myHGG_noDoF.stan") 
  myHGG_noDoF_reparam<-stan_model("/home/aak3/NIST/BTF/BTF_FY21/modelCode/myHGG_noDoF_reparam.stan") 
  
  FirstDatSim=simulateData(N = N,mu = mu,tau = tau,lb = lb,ub=ub)

  saveDat=matrix(c(FirstDatSim$x,FirstDatSim$u),nrow = 1)
  
  if(tau<3 && N<30){
    CIout=oneTrial(dat = FirstDatSim,
                   GGmodel = myHGG_noDoF_reparam,
                   N = N,
                   tau = tau,mu = mu, num.iters = num.iters,alpha = alpha,nu = nu, 
                   howDatSim = howDatSim,rsessionNum = rsessionNum,
                   adaptDeltaVal=theadaptdelta,maxTreeDepthVal=thetreedepth,tdf=thetdf)
  }else{
    CIout=oneTrial(dat = FirstDatSim,
                   GGmodel = myHGG_noDoF,
                   N = N,
                   tau = tau,mu = mu, num.iters = num.iters,alpha = alpha,nu = nu, 
                   howDatSim = howDatSim,rsessionNum = rsessionNum,
                   adaptDeltaVal=theadaptdelta,maxTreeDepthVal=thetreedepth,tdf=thetdf)
    
  }

  CIout$Iteration=1
  CIout$lb=lb
  CIout$ub=ub
  
  for(i in 2:numberOfSims){
    
    newSimDat=simulateData(N = N,mu = mu,tau = tau,lb = lb,ub = ub)
    
    saveDatTemp=matrix(c(newSimDat$x,newSimDat$u),nrow = 1)
    saveDat=rbind(saveDat,saveDatTemp)
    write.csv(saveDat,paste(mypath,"CIouts/",howDatSim,"simulatedData",today,"_tau",tau,"_N",N,"_ub",ub,"_lb",lb,"_",rsession,".csv",sep=""))

    if(tau<3 && N<30){
      newOut=oneTrial(dat = newSimDat,
                     GGmodel = myHGG_noDoF_reparam,
                     N = N,
                     tau = tau,mu = mu, num.iters = num.iters,alpha = alpha,nu = nu, 
                     howDatSim = howDatSim,rsessionNum = rsessionNum,
                     adaptDeltaVal=theadaptdelta,maxTreeDepthVal=thetreedepth,tdf=thetdf)
    }else{
      newOut=oneTrial(dat = newSimDat,
                     GGmodel = myHGG_noDoF,
                     N = N,
                     tau = tau,mu = mu, num.iters = num.iters,alpha = alpha,nu = nu, 
                     howDatSim = howDatSim,rsessionNum = rsessionNum,
                     adaptDeltaVal=theadaptdelta,maxTreeDepthVal=thetreedepth,tdf=thetdf)
      
    }
    
    newOut$Iteration=i
    newOut$lb=lb
    newOut$ub=ub
    CIout=bind_rows(CIout,newOut)
    
    if(is.na(newOut$method[1])){
      print(newOut)
      print(tail(CIout))
    }
    
    if(i %% 10==0) {
      print(paste("iteration",i,sep=" "))
      write.csv(CIout,paste(mypath,"CIouts/",howDatSim,"CIout",today,"_tau",tau,"_N",N,"_ub",ub,"_lb",lb,"_",rsession,".csv",sep=""))
    }
  }
  
  endTime=Sys.time()-startTime
  print(endTime)
  cat(paste("tau is ",tau, ", N is ",N,", endTime is ",endTime,"\n",sep = ""),file = paste(mypath,"CIouts/",howDatSim,"times",today,"_",rsession,".txt",sep=""),append = T)
  
  return("Done!")
  
}

totalTimeStart=Sys.time()

for(paramSetIndex in 1:dim(myparams)[1]){
  oneParamSet=myparams[paramSetIndex,]
  runForManyRepeats(numberOfSims = numberOfSims_set,
                    N = oneParamSet$N_set,
                    tau = oneParamSet$tau_set,
                    mu = mu_set,num.iters = num.iters_set,
                    alpha = alpha_set,nu = nu_set,howDatSim = howDatSim_set,
                    rsessionNum = rsession,
                    theadaptdelta = oneParamSet$adaptDeltaVals,
                    thetreedepth = oneParamSet$maxTreeDepthVals,
                    thetdf = oneParamSet$tdf_set,
                    today = rundate,
                    lb = oneParamSet$lb,
                    ub=oneParamSet$ub)
}

totalTimeEnd=Sys.time()
print(paste("Total time = ",totalTimeEnd-totalTimeStart))
