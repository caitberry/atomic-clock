# source(file = "/home/aak3/NIST/atomic-clock/Code/ComparisonAnalysis2025/DarkUncertaintyMethods/coverageProbExp_withParallel_normal.R")
# .libPaths("/home/aak3/R/x86_64-redhat-linux-gnu-library/3.6/")

rsession=1

library(dplyr)
library(readr)

mypath="/home/aak3/NIST/atomic-clock/Code/ComparisonAnalysis2025/DarkUncertaintyMethods/"

############################################################################
### set parameters based on data
############################################################################

# ratiolab="AlSr"
# ratiolab="AlYb"
ratiolab="YbSr"

ratiodf <- read_csv(paste("/home/aak3/NIST/atomic-clock/Data/ClockComparison2025/BayesianAnalysisData/ErYb_",ratiolab,"_data.csv",sep=""))
# df_YbSr <- read_csv(paste("/home/aak3/NIST/atomic-clock/Data/ClockComparison2025/BayesianAnalysisData/ErYb_YbSr_data.csv")

measurements <- ratiodf$offset
uncertainties <- ratiodf$statistical_unc

N <- length(measurements)

######
mu_set=mean(measurements)

lblabel_set=paste(ratiolab,"_minunc",sep="")
ublabel_set=paste(ratiolab,"_maxunc",sep="")
  
numberOfSims_set=500 ### make bigger

N_set=c(N,N+20,100)
tau_set=c(3,10)

allparams=expand.grid(N_set=N_set,tau_set=tau_set) #want every combo of N and tau

allparams$lb=min(uncertainties)
allparams$ub=max(uncertainties)

myparams=allparams %>%
  mutate(adaptDeltaVals = #########right now, setting these based on N
           ifelse(N_set == 9, .99,
                  ifelse(N_set== 13, .95,
                         ifelse(N_set== 29, .8,
                                ifelse(N_set== 33, .8,
                                       .8))))) %>%
  mutate(maxTreeDepthVals = #########right now, setting these based on N
           ifelse(N_set== 9, 17, 
                  ifelse(N_set== 13, 15,
                         ifelse(N_set== 29, 10,
                                ifelse(N_set== 33, 10,
                                       10))))) %>%
  mutate(tdf_set =
           ifelse(N_set == 5 & tau_set %in%c(1,2), 2,
                  4)) #set prior sd to 4 for all except small N and small tau case, then too big


##############################################################################################
### don't change this (often)
##############################################################################################

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

# # numberOfSims=numberOfSims_set
# # N=N_set[1]
# # tau=tau_set[1]
# # mu=mu_set
# # num.iters=num.iters_set
# # alpha=alpha_set
# # nu=nu_set
# # howDatSim=howDatSim_set
# # # rsession=1
# # # rsessionNum=1
# # theadaptdelta=.99
# # thetreedepth=15
# # thetdf=4
# # today=rundate
# # lb=.1
# # ub=1
# # # 
# 
# paramSetIndex=1
# oneParamSet=myparams[paramSetIndex,]
# numberOfSims = numberOfSims_set
# N = oneParamSet$N_set
# tau = oneParamSet$tau_set
# mu = mu_set
# num.iters = num.iters_set
# alpha = alpha_set
# nu = nu_set
# howDatSim = howDatSim_set
# rsessionNum = rsession
# theadaptdelta = oneParamSet$adaptDeltaVals
# thetreedepth = oneParamSet$maxTreeDepthVals
# thetdf = oneParamSet$tdf_set
# today = rundate
# lb = oneParamSet$lb
# ub=oneParamSet$ub
# 


runForManyRepeats=function(numberOfSims,N,tau,mu,num.iters,alpha,nu,howDatSim,
                           rsessionNum,theadaptdelta,thetreedepth,thetdf,today,lb,ub,lblabel,ublabel){

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
    write.csv(saveDat,paste(mypath,"CIouts/",howDatSim,"simulatedData",today,"_tau",tau,"_N",N,"_ub",ublabel,"_lb",lblabel,"_",rsession,".csv",sep=""))

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
      write.csv(CIout,paste(mypath,"CIouts/",howDatSim,"CIout",today,"_tau",tau,"_N",N,"_ub",ublabel,"_lb",lblabel,"_",rsession,".csv",sep=""))
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
                    ub=oneParamSet$ub,
                    lblabel = lblabel_set,
                    ublabel = ublabel_set)
}

totalTimeEnd=Sys.time()
print(paste("Total time = ",totalTimeEnd-totalTimeStart))
