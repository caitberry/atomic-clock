library(readr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(openxlsx)
library(fields)

library(tidyverse) #data wrangling
library(RSpectra) #eigenvalue solving

mtse=modules::use("/home/aak3/NIST/atomic-clock/Functions.R")


### data analysis steps
# 1. read in data

# read.csv("Data/2024-04-24-ErYb-AlSr.csv",header = F,sep="",na.strings = "NaN")

# need to create dat vector with NAs to send to 2.


########################################################################################
# 2. concatenate small gaps
########################################################################################

smushSmallGaps=function(dat){
  
  missing.indices <- which(is.na(dat))
  differenced <- diff(missing.indices)
  jumps <- which(differenced>1)

  new.dat <- dat
  new.dat[is.na(dat)] <- Inf
  runs <- rle(new.dat)
  NA_period_length <- runs$lengths[runs$values == Inf]
  
  remove.indices <- c()
  start <- 0
  for(i in 1:length(NA_period_length)){
    if(NA_period_length[i]< 100){
      remove.indices <- append(remove.indices,missing.indices[(1 + start):sum(NA_period_length[1:i])])
      start <- sum(NA_period_length[1:i])
    }
    if(NA_period_length[i]>100){
      start <- sum(NA_period_length[1:i])
    }
    
  }
  
  #take out small gaps
  new.dat <- dat[-c(remove.indices)]
  
  return(new.dat)
  
}

#run this to get 
# x.t_all=new.df
# N_tot=length(x.t_all)
# t.vec_all <- 1:N_tot


formatDat=function(dat){
  x.t_all=smushSmallGaps(dat)
  N_tot=length(x.t_all)
  t.vec_all <- 1:N_tot
  
  cleanedDF=data.frame(x.t=x.t_all,t.vec=t.vec_all)
  cleanedDF=cleanedDF%>%mutate(missing=is.na(x.t))
  cleanedDF=filter(cleanedDF,missing==F)
  
  x.t=cleanedDF$x.t
  t.vec <- cleanedDF$t.vec
  N=length(x.t)
  
  return(list(x.t=x.t, t.vec=t.vec, N=N))
}

########################################################################################
# 3. get spectral estimate
# 3a) try a variety of W and K values and check eigenvalues to make sure you've made good selections
# can start with 8/N or 12/N for W
########################################################################################


########################################################################################
# 4. Look at data, spectral estimate, and tapers
########################################################################################

library(RColorBrewer)
# Define the number of colors you want

stuffForPlots=function(myK,x.t,t.vec,specRes){
  nb.cols <- myK+1
  mycolors <- colorRampPalette(brewer.pal(8, "Accent"))(nb.cols)

  resDF=data.frame(n.fourier=length(specRes$freq),
                   freq=specRes$freq,
                   spectrum=specRes$spec.hat)
  
  p1=ggplot(resDF,aes(freq,spectrum,col=factor(interaction(n.fourier))))+
    geom_line()
    # scale_y_log10()+
    # scale_x_log10()
  ## look at tapers
  
  taperDF=data.frame(V.mat$tapers)
  # colnames(taperDF)=c("1","2","3","4","5","6","7","8","9","10","11","12","13",)[1:myK]
  colnames(taperDF)=as.character(1:myK)
  taperDF$t=t.vec
  
  dataDF=data.frame(t=t.vec,
                    value=x.t)
  
  alltData=data.frame(t=min(taperDF$t):max(taperDF$t))
  
  taperDF=merge(alltData,taperDF,by="t",all.x = T)
  
  dataDF=merge(alltData,dataDF,by="t",all.x = T)
  dataDF$type="Data"
  dataDF$Taper="Data"
  
  taperDFlong=melt(taperDF,id.vars = "t",variable.name = "Taper")
  taperDFlong$type="Tapers"
  
  allTaperDat=bind_rows(taperDFlong,dataDF)
  
  p2=ggplot(allTaperDat,aes(t,value,col=Taper))+
    geom_line()+
    facet_wrap(~type,nrow = 2,scales = "free_y",
               strip.position = "left", 
               labeller = as_labeller(c(Data = "Clock Ratio Data", Tapers = "Tapers") ) )+
    # scale_color_brewer(palette = "Accent")+
    scale_color_manual(values=mycolors)+
    # theme(axis.title.y=element_blank())+
    theme(legend.position = "none",legend.title=element_blank(),strip.background = element_blank(),
          strip.placement = "outside")+
    ylab(NULL)
  
  return(list(specPlot=p1,taperPlot=p2))
  
}




########################################################################################
# 5. calculate avar estimate with spectrum for a series of tau values
########################################################################################






########################################################################################
# 6. calculate avar estimate using old method for same tau series
########################################################################################





########################################################################################
# 7. plot both with uncertainties
########################################################################################
