library(tidyr)
library(ggplot2)
library(dplyr)

mtse=modules::use("Functions.R")
source("Functions_SpectrumCovariance.R")

## saving the date to label file outputs
runDate=format(Sys.Date(),"%m%d%y")
numberOfSimulations=300


N <- 1000  # length of the time series
ar_coeffs <- c(0.5, -0.4, 0.3, -0.2)  # AR(4) coefficients

t.vec <- 1:N

# # Simulate AR(4) process
# x.t <- arima.sim(model = list(ar = ar_coeffs), n = N)
# 
# omitted=c(100:300)
# x.t[omitted]=NA
# t.vec[which(is.na(x.t))] <- NA

myK=3
W1factor=5
V.mat3 <- mtse$get_tapers(t.vec, W = W1factor/N, K = myK) 
taperMat=V.mat3$tapers

N.fourier=N/2+1

max.lag.acf=N-1
############################################

MTSEdf=data.frame()
covlist=list()
for(i in 1:numberOfSimulations){
  # print(i)
  # x.t <- rnorm(N)  # Replace with your time series data
  x.t <- arima.sim(model = list(ar = ar_coeffs), n = N)
  # x.t[omitted]=NA
  
  MTSE_full <- mtse$MT_spectralEstimate_fft(x.t, taperMat) 
  MTSE_oneiter_df=data.frame(rep=i,freq=MTSE_full$freqs,spectrum=MTSE_full$spectrum)
  MTSEdf=bind_rows(MTSEdf,MTSE_oneiter_df)
  
  # covlist[[i]]=mtse$spec_cov.mat_slow(x.t = x.t,t.vec = t.vec,taperMat = V.mat3$tapers,N.fourier = N.fourier,isWhite = F,max.lag.acf = max.lag.acf)

}

# average_cov_matrix <- Reduce("+", covlist) / length(covlist)
# 
# matrix_array <- array(unlist(covlist), dim = c(N.fourier, N.fourier, numberOfSimulations))
# variance_matrix <- apply(matrix_array, c(1, 2), sd)


test=mtse$spec_cov.mat_slow(x.t = x.t,t.vec = t.vec,taperMat = V.mat3$tapers,N.fourier = N.fourier,isWhite = F,max.lag.acf = max.lag.acf)
saveRDS(test,file = "Testing/calcCovMatforN1000.rds")
# average_cov_matrix[1:10,1:10]


# MTSEdf %>% group_by(freq) %>% summarise(var=var(spectrum))
diag(test)

# freq1res=MTSEdf %>% filter(freq==0.07) %>% select(spectrum)
# freq2res=MTSEdf %>% filter(freq==0.2) %>% select(spectrum)

myfreqs <- seq(0, 0.5, length.out = N.fourier)
myC=matrix(NA,nrow=N.fourier,ncol=N.fourier)

### Question: is it appropriate to calculate the simulated covariance this way?
for(i in 1:N.fourier){
  for(j in 1:N.fourier){
    freq1res=MTSEdf %>% filter(freq==myfreqs[i]) %>% dplyr::select(spectrum)
    freq2res=MTSEdf %>% filter(freq==myfreqs[j]) %>% dplyr::select(spectrum)
    
    myC[i,j]=cov(freq1res,freq2res)
  }
}

# myC[1:10]
# test[1:10]
# plot(myC[1:51])
# points(Re(test[1:51]),col="blue")
# 
# plot(myC[1:51]-Re(test[1:51]))

howMany=500
myCs=data.frame(index=1:howMany,C=myC[1:howMany],type="From spec reps")
onecalcC=data.frame(index=1:howMany,C=Re(test)[1:howMany],type="one calc C")
# avgCalcC=data.frame(index=1:howMany,C=Re(average_cov_matrix)[1:howMany],type="avg calc C")
# compCs=bind_rows(myCs,onecalcC,avgCalcC)
compCs=bind_rows(myCs,onecalcC)
ggplot(compCs,aes(index,C,color=type))+
  geom_point(alpha=.5)
  # scale_y_log10()
  # facet_wrap(~index,scales = "free")


ggplot(compCs %>% group_by(index) %>% summarise(diff=diff(C)),aes(index,diff))+
  geom_point()



plot(diag(myC)) # calc from sim data
points(diag(Re(test)),col="green") # calculated from out function

plot(diag(Re(test))/diag(myC))
plot(diag(myC)/diag(Re(test)))
1/myK
