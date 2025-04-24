# folderLocation="C:/Users/aak3/Documents/atomic-clock/"
folderLocation="/home/aak3/NIST/atomic-clock/"

source(file = paste(folderLocation,"Functions_SpectrumCovariance.R",sep="")) 
source(file = paste(folderLocation,"Code/analysisFunctions.R",sep=""))

source(file = paste(folderLocation,"Code/ComparisonAnalysis2025/0_data_load.R",sep=""))
source(file = paste(folderLocation,"Code/ComparisonAnalysis2025/1_EDA.R",sep=""))

unique(AlSr_df$date)

# oneRatioDF=YbSr_df
# oneRatioDF=AlYb_df
oneRatioDF=AlSr_df
# for (i in 1:length(unique(oneRatioDF$date))){
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
  
  myDatSpecEst=spectralEstWithUnc(x.t = x.t-mean(x.t),t.vec=t.vec,N.fourier = floor(N/2) + 1, ###############SHOULD I BE REMOVING THE MEAN?????? IT CHANGES THE RESULTS TO LINE UP WITH THE ls RESULTS
                                  ######################################################################## ended the day by rerunning these
                                  
                                  numTapers = myK,calcCov = F,
                                  myW = myW)#,isWhite = T)#,acf.lag = 4,numCores=80)
  myDatSpecEst$e.values
  print(Sys.time()-startTime)
  
  today=format(Sys.Date(),format="%b%d")
  ratdat=format(ratiodate,format="%m%d")
  ratiolab=unique(allDat$ratio)
  # saveRDS(object = myDatSpecEst,file = paste(folderLocation,"Results/ClockComp2025/spectralEstFor",ratiolab,ratdat,"_",today,".Rds",sep=""))

  specRes=myDatSpecEst
  resDF=data.frame(freq=specRes$freq,
                   spectrum=specRes$spec.hat,
                   method = "MTSE")
  
  
  
  p1=ggplot(resDF,aes(freq,spectrum,col=method))+
    geom_line() +
    scale_y_log10()+
    scale_x_log10()
  p1
  
  
  ########################################################################################
  # 3b. get LS spectral estimate
  ########################################################################################
  n.test=t.vec[length(t.vec)]-t.vec[1]
  N.fourier = floor(n.test/2) + 1
  myfreqs <- seq(0,0.5, length.out = N.fourier)
  length(myfreqs)
  
  LSfreq=myfreqs#seq(0,.5,length.out=100)
  # t.vec_shift <-  (t.vec - min(t.vec)) / max(t.vec)  # scale between 0 and 1
  LSperDF=data.frame(freq=LSfreq,#specRes$freq,
                     spectrum=mtse$lomb_scargle_processedDat(x.t = x.t,t.vec = t.vec,f = LSfreq),
                     method="LS_mine")
  # # LSfreq=seq(0,.5,length.out=100)
  # # LSperDF=data.frame(freq=LSfreq,#specRes$freq,
  # #                    spectrum=mtse$lomb_scargle(x.t = dat,f = LSfreq),
  # #                    method="LS")
  # 
  # library(lomb)
  # lombres=lsp(x = x.t,times = t.vec)
  # lombrespress=lsp(x = x.t,times = t.vec,normalize = "press")
  # lombresDF=data.frame(freq=lombrespress$scanned,#specRes$freq,
  #                    spectrum=lombrespress$power*var(x.t),
  #                    method="LS_lombpackage")   
  # # my/cait code is very very close (not exactly same) to results from lomb package if i use normalize = "press" and resulting power*var(x.t)
  # ggplot(lombresDF,aes(freq,spectrum,col=method))+
  #   geom_line() +
  #   scale_y_log10()+
  #   scale_x_log10()
  # # 
  
  lomb_power <- function(t, y, freqs) {
    n <- length(freqs)
    power <- numeric(n)
    
    for (i in seq_along(freqs)) {
      omega <- 2 * pi * freqs[i]
      
      cos2wti <- cos(2 * omega * t)
      sin2wti <- sin(2 * omega * t)
      
      tau <- atan2(sum(sin2wti), sum(cos2wti)) / (2 * omega)
      
      wt <- omega * (t - tau)
      cos_wt <- cos(wt)
      sin_wt <- sin(wt)
      
      C <- sum(y * cos_wt)
      S <- sum(y * sin_wt)
      
      CC <- sum(cos_wt^2)
      SS <- sum(sin_wt^2)
      
      power[i] <- 0.5 * ((C^2 / CC) + (S^2 / SS))
    }
    
    return(power)
  }
  
  power <- lomb_power(t.vec, x.t-mean(x.t), myfreqs)
  
  GPTresDF=data.frame(freq=myfreqs,#specRes$freq,
                       spectrum=power,
                       method="GPTlomb")
  # 
  # # 
  ggplot(LSperDF,aes(freq,spectrum,col=method))+
    geom_line() +
    scale_y_log10()+
    scale_x_log10()
  # 
  bothspec=bind_rows(LSperDF,GPTresDF,resDF)
  
  ggplot(bothspec,aes(freq,spectrum,col=method))+
    geom_line(alpha=.4) +
    scale_y_log10()+
    scale_x_log10()
    # coord_cartesian(xlim = c(0.0001,.001))
  
  head(filter(bothspec,freq>1.299376e-04, freq<.1,method=="LS_mine"))#$spectrum
  head(filter(bothspec,freq>1.299376e-04, freq<.1,method=="LS_lombpackage"))#$spectrum
  head(filter(bothspec,freq>1.299376e-04, freq<.1,method=="GPTlomb"))#$spectrum
  # my answer and GPT answer look the same to me exaclty, package answer is close
  # GPT and mine only matched after i demeaned the data. 
  # do i need to demean mtse data?
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
  
  oldAVARestUncertainty=mtse$avar_CI(CI.level = .68,
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
  
# }


