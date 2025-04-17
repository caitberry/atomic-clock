# folderLocation="C:/Users/aak3/Documents/atomic-clock/"
folderLocation="/home/aak3/NIST/atomic-clock/"

source(file = paste(folderLocation,"Functions_SpectrumCovariance.R",sep="")) 
source(file = paste(folderLocation,"Code/analysisFunctions.R",sep=""))

source(file = paste(folderLocation,"Code/ComparisonAnalysis2025/0_data_load.R",sep=""))
source(file = paste(folderLocation,"Code/ComparisonAnalysis2025/1_EDA.R",sep=""))

unique(AlSr_df$date)

# oneRatioDF=YbSr_df
# oneRatioDF=AlYb_df
# oneRatioDF=AlSr_df

allRatios=bind_rows(AlSr_df,AlYb_df,YbSr_df)

for(j in 1:3){
  oneRatioDF=filter(allRatios, ratio==unique(allRatios$ratio)[j])
  for (i in 1:length(unique(oneRatioDF$date))){
      
    ratiodate=unique(oneRatioDF$date)[i]
    print(ratiodate)
    allDat=filter(oneRatioDF, date == ratiodate)
    
    ########################################################################################
    # 1. read in data
    ########################################################################################
    dat=allDat$V2
    
    ########################################################################################
    # 2. concatenate small gaps
    ########################################################################################
    formattedData=formatDat(dat)
    
    x.t=formattedData$x.t
    t.vec=formattedData$t.vec
    N=formattedData$N
    
    ########################################################################################
    # 3. get spectral estimate
    # 3a) try a variety of W and K values and check eigenvalues to make sure you've made good selections
    # can start with 8/N or 12/N for W
    ########################################################################################
    myW=8/N
    myK=5
    # V.mat <- mtse$get_tapers(t.vec, W = myW, K = myK) 
    # print(V.mat$e.values)
    startTime=Sys.time()
    
    myDatSpecEst=spectralEstWithUnc(x.t = x.t-mean(x.t),t.vec=t.vec,N.fourier = floor(N/2) + 1,
                                    numTapers = myK,calcCov = F,
                                    myW = myW)#,isWhite = T)#,acf.lag = 4,numCores=80)
    print(Sys.time()-startTime)
    
    today=format(Sys.Date(),format="%b%d")
    ratdat=format(ratiodate,format="%m%d")
    ratiolab=unique(allDat$ratio)
    saveRDS(object = myDatSpecEst,file = paste(folderLocation,"Results/ClockComp2025/spectralEstFor",ratiolab,ratdat,"_",today,".Rds",sep=""))
    
    ########################################################################################
    # 5. calculate avar estimate with spectrum for a series of tau values
    ########################################################################################
    taus <- 2^(0:9)
    taus <- taus[taus<floor(N/3)]
    
    specAVARest=mtse$AVAR_spec(spectral_est = myDatSpecEst$spec.hat,taus = taus,calcUnc = F)#,Cov.mat = myDatSpecEst$Cov.mat)
    
    ########################################################################################
    # 6. calculate avar estimate using old method for same tau series
    ########################################################################################
    
    oldAVARest=mtse$overlapping_avar_fn(y = x.t,m = taus)
    
    oldAVARestUncertainty=mtse$avar_CI(CI.level = .95,
                                       noise_type = "white noise", 
                                       avar_type = "chisquared", 
                                       avars = oldAVARest, 
                                       taus=taus,
                                       N=N)
    ########################################################################################
    # 7. plot both with uncertainties
    ########################################################################################
    
    avarDF=data.frame(avar=c(oldAVARest,specAVARest$avar),
                      avarLower=c(oldAVARestUncertainty$lower,specAVARest$avar-sqrt(specAVARest$avarVar)),
                      avarUpper=c(oldAVARestUncertainty$upper,specAVARest$avar+sqrt(specAVARest$avarVar)),
                      tau=rep(taus,2),
                      Method=rep(c("Current","Spectral"),each=length(taus)),
                      Data = paste(ratiolab,ratdat,sep=""))
    
    saveRDS(object = avarDF,file = paste(folderLocation,"Results/ClockComp2025/avarDFfor",ratiolab,ratdat,"_",today,".Rds",sep=""))
    
  }
  
}

# colnames(allDat)=c("MJD","AlSr_20250227")
# ggplot(allDat,aes(MJD,AlSr_20250227))+
#   geom_point()
# ggplot(filter(allDat,MJD<60733.85 & MJD>60733.84),
#        aes(MJD,AlSr_20250227))+
#   geom_point()

########################################################################################
# 4. Look at data, spectral estimate, and tapers
########################################################################################
plotOutput=stuffForPlots(myK,x.t,t.vec,specRes=myDatSpecEst)
plotOutput$specPlot
plotOutput$taperPlot



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
