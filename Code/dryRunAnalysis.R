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



###############################################################################################
###############################################################################################

# 1. read in data
dat2024_04_24_ErYb_AlSr <- read.csv("Data/clockRatioTimeseries_first2024dryRun/2024-04-24-ErYb-AlSr.csv",header = F,sep="",na.strings = "NaN")
dim(dat2024_04_24_ErYb_AlSr)
colnames(dat2024_04_24_ErYb_AlSr)=c("MJD","FracDiff")

plot(dat2024_04_24_ErYb_AlSr$MJD,dat2024_04_24_ErYb_AlSr$FracDiff)
dataInDF=data.frame(MJD=dat2024_04_24_ErYb_AlSr$MJD,y=dat2024_04_24_ErYb_AlSr$FracDiff-mean(dat2024_04_24_ErYb_AlSr$FracDiff,na.rm = T))

ggplot(dataInDF,aes(MJD,y))+
  geom_point()

ggplot(dat2024_04_24_ErYb_AlSr,aes(MJD,FracDiff))+
  geom_point()

dat2024_04_24_ErYb_AlSr$seconds=1:dim(dat2024_04_24_ErYb_AlSr)[1]


dat=dat2024_04_24_ErYb_AlSr%>%mutate(missing=is.na(FracDiff),
                                     date=convertToDateTime(MJD, origin = "1858-11-17",tz="MDT"))

plot(dat$MJD,dat$FracDiff)

# dat=filter(dat, missing==F)
# t.vec <- dat$seconds
# x.t=dat$FracDiff#-mean(dat$FracDiff)
# N=length(x.t)

# 2. concatenate small gaps
## need to remove the tiny gaps as they mess up the tapers, this is an assumption we are making 

####EDA on missing points

sum(is.na(dat))
df=dat$FracDiff
missing.indices <- which(is.na(df))
missing.indices[1:20]
differenced <- diff(missing.indices)
differenced[1:20]
sum(differenced > 1)
jumps <- which(differenced>1)
differenced[1220:1230]

missing.indices[1:10]
differenced[1:10]

new.df <- df
new.df[is.na(df)] <- Inf
runs <- rle(new.df)
runs
NA_period_length <- runs$lengths[runs$values == Inf]

NA_period_length


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
new.df <- df[-c(remove.indices)]

#look at gaps now
new.df2 <- new.df
new.df2[is.na(new.df)] <- Inf
runs2 <- rle(new.df2)
runs2
NA_period_length2 <- runs2$lengths[runs2$values == Inf]

NA_period_length2


### define the data


x.t_all=new.df
# # ##############################testing other x.t
# # x.t_all=rnorm(1000)
# # x.t_all=((new.df-mean(new.df,na.rm = T))/sd(new.df,na.rm = T))[1:1000]
# x.t_all=new.df*10^16
N_tot=length(x.t_all)
t.vec_all <- 1:N_tot

plot(t.vec_all,x.t_all)


cleanedDF=data.frame(x.t=x.t_all,t.vec=t.vec_all)
cleanedDF=cleanedDF%>%mutate(missing=is.na(x.t))
cleanedDF=filter(cleanedDF,missing==F)

x.t=cleanedDF$x.t
t.vec <- cleanedDF$t.vec
N=length(x.t)

plot(t.vec,x.t)

length(x.t)
length(t.vec)

# 3. get spectral estimate
# 3a) try a variety of W and K values and check eigenvalues to make sure you've made good selections
# can start with 8/N or 12/N for W
V.mat <- mtse$get_tapers(t.vec, W = 4/N*3, K = 10) 
V.mat$e.values
test=mtse$spectralEstWithUnc(x.t = x.t,t.vec=t.vec,N.fourier = floor(N/2) + 1,#100,
                        numTapers = 10,calcCov = F,
                        myW = 4/N*3)
test$V.mat$e.values ####NEED TO ADD THIS ABILITY TO THE FUNCTION
dim(test$Cov.mat)

# # save test
# today=format(Sys.Date(),format="%b%d")
# resName=paste("2024-04-24-ErYb-AlSr",today,"firstTry",sep="_")
# 
# saveRDS(test,file = paste("/home/aak3/NIST/atomic-clock/Results/resultsFor",resName,".Rds",sep=""))

# 4. Look at data, spectral estimate, and tapers
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 10+1
mycolors <- colorRampPalette(brewer.pal(8, "Accent"))(nb.cols)

datDF=data.frame(x.t=x.t-mean(x.t),t.vec=t.vec)
ggplot(datDF,aes(t.vec,x.t))+
  geom_point()
  
resDF=data.frame(ratio="test",
                n.fourier=length(test$freq),
                freq=test$freq,
                spectrum=test$spec.hat)

ggplot(resDF,aes(freq,spectrum,col=factor(interaction(n.fourier))))+
  geom_line()+
  scale_y_log10()+
  scale_x_log10()

## look at tapers

taperDF=data.frame(V.mat$tapers)
colnames(taperDF)=c("1","2","3","4","5","6","7","8","9","10")
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

ggplot(allTaperDat,aes(t,value,col=Taper))+
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


# 5. calculate avar estimate with spectrum for a series of tau values


taus <- 2^(0:9)
taus <- taus[taus<floor(N/3)]
  
specAVARest=mtse$AVAR_spec(spectral_est = test$spec.hat,taus = taus,calcUnc = F,Cov.mat = test$Cov.mat)

# 6. calculate avar estimate using old method for same tau series

plot(x.t)

oldAVARest=mtse$overlapping_avar_fn(y = x.t,m = taus)

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

