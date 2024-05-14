# Title: Functions Module
# File name: Functions.R
# Description: This file contains all of the functions relevant to 
#   calculating the spectral estimate and allan variance of time series data
#   possibly with gaps

require(modules)

# libraries ###########################################################################
import("tidyverse", #data wrangling
       "RSpectra", #eigenvalue solving
       "fields", #dist.mat
       "dplyr") #bind_rows

# export("get_tapers", "MT_spectralEstimate") #only use if we don't want to export every function

# Functions ###########################################################################

## Calculating Missing Data Tapers ####################################################
#inputs: (->) t.n = time points of length N (possibly with NA values), 
#        (->) W = analysis half bandwidth, 
#        (->) K = number of tapers
#output: (<-) L x K matrix of tapers where 
#               L = length of time series without missing values, 
#               K = number of tapers

get_tapers <- function(t.n, W, K){
  dist.mat <- rdist(na.omit(t.n))
  
  #create the A' matrix (Chave 2019 equation (22))
  A.prime <- (1/(pi*dist.mat))*sin(2*pi*W*dist.mat)
  A.prime[row(A.prime) == col(A.prime)] <- W*2
  print("A matrix computed")
  
  eigdec <- eigs(A.prime, k = K, which = "LM")
  eig_vecs <- eigdec$vectors #get only the vectors
  print("tapers computed")
  
  return(list("tapers" = eig_vecs, "e.values" = eigdec$values))
} 


## Calculating Missing Data MTSE #####################################################

#inputs:  (->) X.t = time series of length N with any missing values 
#                   and length L without, 
#         (->) V.mat = L X K dimension taper matrix
#outputs: (<-) freqs = fourier frequencies
#         (<-) spectrum = spectral estimate

MT_spectralEstimate <- function(X.t, V.mat){
  im <- complex(real = 0, imaginary = 1)
  X.t <- X.t - mean(X.t, na.rm = TRUE) #demean
  N.long <- length(X.t)
  t.n <- 1:N.long
  missing.indices <- which(is.na(X.t))
  t.n[which(is.na(X.t))] <- NA
  t.n_m <- rdist(na.omit(t.n))[1,]
  
  ##use tapers to generate spectral estimate
  N <- length(na.omit(t.n))
  S.x.hat <- rep(NA, times = floor(N/2) + 1)
  freqs <- seq(0,0.5, length.out = floor(N/2) + 1)
  K <- dim(V.mat)[2]
  
  for(j in 1:length(freqs)){
    k.vec <- rep(NA,times = K)
    for(k in 1:K){
      inner.sum <- sum(V.mat[,k]*na.exclude(X.t)*exp(-im*2*pi*freqs[j]*t.n_m))
      k.vec[k] <- abs(inner.sum)^2
    }
    S.x.hat[j] <- mean(k.vec)
  }
  return(list("spectrum" = S.x.hat, "freqs" = freqs))
}



## Calculating MTSE with FFT ##########################################################

#inputs:  (->) X.t = time series of length N with any missing values 
#                   and length L without, 
#         (->) V.mat = L X K dimension taper matrix
#outputs: (<-) freqs = fourier frequencies
#         (<-) spectrum = spectral estimate


MT_spectralEstimate_fft <- function(X.t, V.mat){
  
  ##use tapers to generate spectral estimate
  N <- length(na.exclude(X.t))
  N.fourier <- floor(N/2) + 1
  S.x.hat <- rep(NA, times = N.fourier)
  freqs <- seq(0,0.5, length.out = N.fourier)
  K <- dim(V.mat)[2]
  S.k.mat <- matrix(NA,nrow = K, ncol = N.fourier)
  
  for(k in 1:K){
    spec.vec <- fft(taperMatrix[,k]*na.exclude(X.t))[1:N.fourier]
    S.k.mat[k,] <- abs(spec.vec)^2
  }
  
  S.x.hat <- apply(S.k.mat, MARGIN = 2, FUN = mean)
  
  return(list("spectrum" = S.x.hat, "freqs" = freqs))
}

## AVAR Calculation (Regular) #########################################################

#inputs:  (->) y = time series of length n without any missing values
#         (->) tau = averaging time        
#output:  (<-) AVAR estimate for given tau


avar_fn <- function(y,tau){
  n=length(y)
  
  div=seq(1,n,by = tau)
  
  M=length(div)-1 #number of groups
  
  groupmeans = numeric(M)
  for(i in 1:M){
    groupmeans[i]=mean(y[div[i]:(div[i+1]-1)])
  }
  
  1/(2*(M-1)) * sum(diff(groupmeans)^2)
}

## OAVAR Calculation (Overlapping) ####################################################

#inputs:  (->) y = time series of length n without any missing values
#         (->) m = averaging time        
#output:  (<-) Overlapping AVAR estimate for given m

overlapping_avar_fn <- function(y,m){
  
  M=length(y)
  
  outer.sum = 0
  for(j in 1:(M-2*m+1)){
    sum = 0
    for(i in j:(j + m - 1)){
      sum = sum+ y[i+m] - y[i]
    }
    outer.sum = outer.sum + sum^2
  }
  
  out <- 1/(2*m^2*(M - 2*m + 1))*outer.sum
  
  return(out)
}

## True AVAR line for ARFIMA(0,d,0) ###################################################

#inputs:  (->) N.tau = max tau value you would like calculated
#         (->) d = parameter of ARFIMA process
#         (->) sig.2.a = variance of white noise process in ARFIMA process
#output:  (<-) vector of length N.tau-1 with  

tavar_ARFIMA <- function(N.tau,d, sig.2.a){
  rho.vec <- tacvfARFIMA(phi = 0, theta = 0, dfrac = d, maxlag = 2*N.tau)
  corr.vec <- rho.vec/max(rho.vec) #normalize
  taus <- 2:N.tau
  
  sum.vec <- rep(NA, times = N.tau-1)
  for(k in 2:N.tau){
    #print(k)
    total <- 0
    for(i in 1:(k - 1)){
      #print(i) 
      total = total + i*(2*corr.vec[k-i + 1] - corr.vec[i + 1] - corr.vec[2*k-i + 1])
    }
    sum.vec[k-1] <- total
  }
  numerator <-(taus*(rep(corr.vec[1],times = N.tau-1) - corr.vec[3:(N.tau+1)]) + sum.vec)*gamma(1-2*d)
  denom <- (taus*gamma(1-d))^2
  return(numerator/denom)
}


## Calculating the Lomb Scargle Periodogram ###########################################

#inputs:  (->) x.t = vector of data (possibly with NA values)
#         (->) f = frequencies at which we calculate the lomb-scargle periodogram
#output:  (<-) spectral estimate


#lomb-scargle function
lomb_scargle <- function(x.t,f){
  
  ## calculates the Lomb-Scargle Periodogram for data x.t at frequencies f##
  
  N <- length(x.t)
  L <- length(f)
  t.vec <- 1:N
  t.vec[which(is.na(x.t))] <- NA
  x.missing <- na.exclude(x.t)
  t.missing <- na.exclude(t.vec)
  
  lsperio <- rep(NA, times = L)
  
  for(i in 1:L){
    x.centered <- x.missing - mean(x.missing)
    x.var <- var(x.missing)
    tau.value <- tau.shift(f[i], t.missing)
    c.vec <- cos(2*pi*(f[i]*(t.missing - tau.value)))
    s.vec <- sin(2*pi*(f[i]*(t.missing - tau.value)))
    lsperio[i] <- (1/(2*x.var))*((x.centered%*%c.vec)^2/sum(c.vec^2) + 
                                   (x.centered%*%s.vec)^2/sum(s.vec^2))
  }
  
  return(lsperio)
}


#tau function (local function needed for lomb_scargle())

#inputs:  (->) f = frequencies at which we calculate the lomb-scargle periodogram
#         (->) t = 
#output:  (<-) 

tau.shift <- function(f,t){
  (1/(4*pi*f))*atan(sum(sin(4*pi*f*t))/sum(cos(4*pi*f*t)))
}

## Uncertainty: Spectral-based Method ################################################

## Spectral Estimate and Uncertainty ################################################
#inputs:  (->) x.t = vector of data (possibly with NA values)
#         (->) numTapers = number of tapers
#output:  (<-) spectral estimate with covariance matrix

spectralEstWithUnc=function(x.t,numTapers){
  N <- length(x.t)
  t.vec <- 1:N
  t.vec[which(is.na(x.t))] <- NA
  
  V.mat <- get_tapers(t.vec, W = 7/N, K = numTapers)
  MTSE_full <- MT_spectralEstimate(x.t, V.mat$tapers) #new function, calculates just spectrum
  N <- length(t.vec)
  N.fourier <- floor(N/2) + 1
  freq <- seq(0,0.5, length.out = N.fourier)
  
  delta.f <- freq[2] #interval spacing between frequencies, needed for spectral avar calculation
  
  ### calculate the covariance matrix 
  Cov.mat_chave <- matrix(NA, nrow = N.fourier, ncol = N.fourier)
  
  for(i in 1:N.fourier){
    j = 1
    while(j <= i){
      Cov.mat_chave[i,j] <- norm(Conj(t(V.mat$tapers*exp(-im*2*pi*freq[i]*t.vec)*(1/sqrt(numTapers))))%*%(V.mat$tapers*exp(-im*2*pi*freq[j]*t.vec)*(1/sqrt(numTapers))), type = "2") 
      j = j+1
    }
  }
  
  Cov.mat_chave[upper.tri(Cov.mat_chave)] <- t(Cov.mat_chave)[upper.tri(Cov.mat_chave)]
  
  return(list(freq=freq,
              spec.hat=MTSE_full$spectrum,Cov.mat=Cov.mat_chave))
  
}

## G_tau(f) transfer function #########################################################
#inputs:  (->) f = frequencies at which we calculate the transfer function
#         (->) tau = single tau value to calculate transfer function for
#output:  (<-) vector giving value of transfer function for tau and all frequencies

transfer.func <- function(f,tau){
  4*sin(pi*f*tau)^4/(tau*sin(pi*f))^2
}

## AVAR Calculation (Spectral) #########################################################
# (eq'n at bottom of p. 70 of reappraisal) 
#inputs:  (->) spectral_est = spectral estimate (as a vector)
#         (->) taus = taus (as a vector) where you want the AVAR calculated
#         (->) Cov.mat_chave = covariance matrix estimate. Provide to get uncertainty estimate for avar
#output:  (<-) spectral AVAR estimate with variance estimate (if covariance of the spectrum is provided)

AVAR_trfunc <- function(spectral_est, taus,Cov.mat_chave=NA){
  f <- seq(0,0.5,length.out = length(spectral_est))
  
  cov.mat=avar=numeric(length(taus))
  
  for(i in 1:length(taus)){
    G.vec <- transfer.func(f,tau = taus[i]) 
    G.vec[1] <- 0 
    
    avar[i]=f[2]*sum(G.vec*spectral_est)
    
    if(!is.na(Cov.mat_chave)){
      #calculate variance for the AVAR estimate at the given tau
      cov.mat[i] <- t(G.vec)%*%(Cov.mat_chave)%*%G.vec*(f[2])^2
    }
    
  }
  return(list(avar=avar,avarVar=cov.mat))
}


## Uncertainty: Overlapping AVAR #####################################################
#inputs: (->) CI.level = desired confidence level for interval 
#        (->) noise_type = assumed noise type for the data. Right now only works for "white noise"
#        (->) avar_type = either 'simple' or 'chisquared', corresponding to intervals described on pages 37/38 of THOFSA
#        (->) avars = allan variance estimates
#        (->) taus = taus (as a vector) where you want the AVAR calculated
#        (->) N = length of the data (with no gaps)
#output: (<-) CI.limits = a tibble with 3 columns, one for tau, one for the lower bound 
#                         and one for the upper bound of the CI
#               
#               

avar_CI <- function(CI.level,noise_type = "white noise", avar_type, avars, taus,N){
  
  a <- (1-CI.level)/2
  s.2=avars
  
  edf <- rep(NA, times = length(taus))
  i=1
  
  if(noise_type == "white noise"){
    for(m in taus){
      edf[i] <- ((3*(N-1)/(2*m)) - (2*(N-2)/N))*(4*m^2)/(4*m^2 + 5)
      i=i+1
    }
  }
  
  if(avar_type == "simple"){
    CI.limits <- bind_rows("lower" = s.2 - s.2/N, "upper" = s.2 +  s.2/N)
  }else if(avar_type == "chisquared"){
    CI.limits <- bind_rows("lower" = s.2*edf/qchisq(1-a,edf),"upper" = s.2*edf/qchisq(a, edf) )
  }else{
    warning("Invalid avar_type, should be 'simple' or 'chisquared'")
  }
  
  CI.limits = bind_cols(data.frame(tau=taus),CI.limits)
  
  return(CI.limits)
}

