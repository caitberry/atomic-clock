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
       "dplyr", #bind_rows
       "stats", #na.omit, na.exlclude
       "arfima",#tacvfARFIMA
       "parallel",
       "Rcpp",
       "RcppArmadillo") 

export("get_tapers", 
       "MT_spectralEstimate_fft",
       "avar_fn_vec", "overlapping_avar_fn_vec", "tavar_ARFIMA",
       "lomb_scargle", "transfer.func",
       "AVAR_spec","avar_CI") #functions to export


# Helper Functions ####################################################################

## Inserting NA values into a vector at specified indices 
#inputs: (->) original_vector = L x 1 vector of non-NA values
#        (->) template_vector = N x 1 vector of values and NAs in the 
#                               correct index placement which we want to match in the output
#output: (<-) vector of size N x 1

.fill_vector_with_na_sprinkled <- function(original_vector, template_vector) {
  # Length of the original vector
  N <- length(original_vector)
  
  # Check that the template vector has the appropriate length
  if (sum(is.na(template_vector)) != (length(template_vector) - N)) {
    stop("The number of NA values in the template vector must be equal to the difference in lengths.")
  }
  
  # Initialize a counter for the original vector
  original_index <- 1
  
  # Fill in the template vector
  for (i in 1:length(template_vector)) {
    if (is.na(template_vector[i])) {
      next
    } else {
      template_vector[i] <- original_vector[original_index]
      original_index <- original_index + 1
    }
  }
  
  return(template_vector)
}


# Main Functions ###########################################################################

## Calculating Missing Data Tapers ####################################################
#inputs: (->) t.n = time points of length N (possibly with NA values), 
#        (->) W = analysis half bandwidth, 
#        (->) K = number of tapers
#output: (<-) N x K matrix of tapers where 
#               N = length of time series with missing values, 
#               K = number of tapers

get_tapers <- function(t.n, W, K){
  dist.mat <- fields::rdist(stats::na.omit(t.n))
  
  #create the A' matrix (Chave 2019 equation (22))
  A.prime <- (1/(pi*dist.mat))*sin(2*pi*W*dist.mat)
  A.prime[row(A.prime) == col(A.prime)] <- W*2
  print("A matrix computed")
  
  eigdec <- RSpectra::eigs(A.prime, k = K, which = "LM")
  eig_vecs <- eigdec$vectors #get only the vectors
  print("tapers computed")
  
  taperMatrixNA <- apply(eig_vecs, MARGIN = 2, FUN = .fill_vector_with_na_sprinkled, template_vector = t.n)
  
  return(list("tapers" = taperMatrixNA, "e.values" = eigdec$values))
} 



## Calculating Missing Data MTSE #####################################################
#inputs:  (->) X.t = time series of length N including any missing values 
#         (->) V.mat = N X K dimension taper matrix
#outputs: (<-) freqs = fourier frequencies
#         (<-) spectrum = spectral estimate



MT_spectralEstimate_fft <- function(X.t, V.mat){
  ##use tapers to generate spectral estimate
  N <- length(X.t)
  X.t[is.na(X.t)] <- 0 #set NA values in X.t to 0
  V.mat[is.na(V.mat)] <- 0 #set NA values in V.mat to 0
  N.fourier <- floor(N/2) + 1
  S.x.hat <- rep(NA, times = N.fourier)
  freqs <- seq(0,0.5, length.out = N.fourier)
  K <- dim(V.mat)[2]
  S.k.mat <- matrix(NA,nrow = K, ncol = N.fourier)
  
  for(k in 1:K){
    spec.vec <- stats::fft(V.mat[,k]*X.t)[1:N.fourier]
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



## AVAR Calculation (Regular, Vectorized) #########################################################

#inputs:  (->) y = time series of length n without any missing values
#         (->) tau = vector of averaging times        
#output:  (<-) AVAR estimate for given tau


avar_fn_vec <- Vectorize(FUN = avar_fn, vectorize.args = "tau")




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


## OAVAR Calculation (Overlapping, Vectorized) ####################################################

#inputs:  (->) y = time series of length n without any missing values
#         (->) m = vector of averaging times        
#output:  (<-) Overlapping AVAR estimate for given m's

overlapping_avar_fn_vec <- Vectorize(FUN = overlapping_avar_fn, vectorize.args = "m")




## True AVAR line for ARFIMA(0,d,0) ###################################################

#inputs:  (->) N.tau = max tau value you would like calculated
#         (->) d = parameter of ARFIMA process
#         (->) sig.2.a = variance of white noise process in ARFIMA process
#output:  (<-) vector of length N.tau-1 with  

tavar_ARFIMA <- function(N.tau,d, sig.2.a){
  rho.vec <- arfima::tacvfARFIMA(phi = 0, theta = 0, dfrac = d, maxlag = 2*N.tau)
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
  x.missing <- stats::na.exclude(x.t)
  t.missing <- stats::na.exclude(t.vec)
  
  lsperio <- rep(NA, times = L)
  
  for(i in 1:L){
    x.centered <- x.missing - mean(x.missing)
    x.var <- stats::var(x.missing)
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



## Spectrum Covariance Matrix Function ######################################################
# Calculates Cov(S.hat(f_i), S.hat(f_j))
#For data of length N, possibly with gaps
#length N.omitted without gaps
#N.fourier is the number of fourier frequencies
#inputs: (->) X.t = Nx1 vector of data, possibly with gaps
#        (->) N.fourier
#        (->) taperMat = N.omitted x K matrix of data tapers calculated from get_tapers() function, 
#                       where K is the number of tapers 
#        (->) isWhite = boolean which is set to true if assuming white noise and false if not
#        (->) acf.lag
#output: (<-) N.fourier x N.fourier matrix where the ijth entry is Cov(S.hat(f_i), S.hat(f_j))





## G_tau(f) transfer function ########################################################
#inputs:  (->) f = frequencies at which we calculate the transfer function
#         (->) tau = single tau value to calculate transfer function for
#output:  (<-) vector giving value of transfer function for tau and all frequencies

transfer.func <- function(f,tau){
  4*sin(pi*f*tau)^4/(tau*sin(pi*f))^2
}


## Spectral-Based AVAR Calculation ##################################################
# (eq'n at bottom of p. 70 of reappraisal) 
#inputs:  (->) spectral_est = spectral estimate (as a vector)
#         (->) taus = taus (as a vector) where you want the AVAR calculated
#         (->) calcUnc = Boolean where T = you want the uncertainty estimate and F = you do not want the uncertainty
#         (->) Cov.mat_chave = covariance matrix estimate. Provide to get uncertainty estimate for avar
#output:  (<-) spectral AVAR estimate with variance estimate (if covariance of the spectrum is provided)

AVAR_spec <- function(spectral_est, taus, calcUnc = F, Cov.mat =NA){
  f <- seq(0,0.5,length.out = length(spectral_est))
  
  cov.mat=avar=numeric(length(taus))
  
  for(i in 1:length(taus)){
    G.vec <- transfer.func(f,tau = taus[i]) 
    G.vec[1] <- 0 
    
    avar[i]=f[2]*sum(G.vec*spectral_est)
    
    if(calcUnc){
      #calculate variance for the AVAR estimate at the given tau
      cov.mat[i] <- t(G.vec)%*%(Cov.mat)%*%G.vec*(f[2])^2
    }
    
  }
  return(list(avar=avar,avarVar=cov.mat))
}

## Spec-based AVAR Variance Calculation #############################################
#inputs:  (->) taus = vector of tau values at which you'd like the uncertainty calculated
#         (->) C.mat = Covariance matrix object from function spec_cov.mat(), an N.fourier x N.fourier sized matrix
#                       where N.fourier is the number of fourier frequencies
#output:  (<-) avar_var = vector of Var(sigma^2(taus)) values to be used to calculate the bars   

AVAR_spec_var <- function(taus, C.mat){
  variance_vec <- rep(NA, times = length(taus))
  freqs <- seq(0,0.5, length.out = dim(C.mat)[1])
  
  for(i in 1:length(taus)){
  G.vec <- transfer.func(taus[i], freqs)
  variance_vec[i] <- t(G.vec)%*%(C.mat)%*%G.vec*(freqs[2])^2
  }
}
  
## Spec-based AVAR Uncertainty Calculation #############################################
#inputs:  (->) CI.level = desired confidence level for interval 
#         (->) taus = vector of tau values at which you'd like the uncertainty calculated
#         (->) avar = vector of sigma^2(taus) values to be used to calculate the bars 
#         (->) avar_var = vector of Var(sigma^2(taus)) values to be used to calculate the bars 
#output: (<-) CI.limits = a tibble with 3 columns, one for tau, one for the lower bound 
#                         and one for the upper bound of the CI on the spectral based AVAR estimate

AVAR_spec_CI <- function(CI.level, taus, avar, avar_var){

  a <- (1-CI.level)/2
  
  CI.limits <- dplyr::bind_rows("lower" = avar + qnorm(a)*sqrt(avar_var),
                                "upper" = avar + qnorm(1-a)*sqrt(avar_var))
  
  CI.limits = dplyr::bind_cols(data.frame(tau=taus),CI.limits)
  
  return(CI.limits)
  
}
  
  
  
# ## Spectral Estimate and Uncertainty ################################################
# #inputs:  (->) x.t = vector of data (possibly with NA values, MT_spectralEstimate_fft excludes them)
# #         (->) t.vec = vector of times (possibly with NA values, but doesn't need them for get_tapers)
# #         (->) N.fourier
# #         (->) numTapers = number of tapers
# #         (->) calcCov = True if you want to calculate the covariance matrix
# #         (->) myW = Analysis half-bandwidth
# #         (->) isWhite
# #         (->) acf.lag
# #output:  (<-) spectral estimate with covariance matrix
# 
# spectralEstWithUnc <- function(x.t,t.vec,N.fourier,numTapers,calcCov=T,myW,isWhite = TRUE,acf.lag=4,numCores){
#   freq <- seq(0,0.5, length.out = N.fourier)
#   
#   delta.f <- freq[2] #interval spacing between frequencies, needed for spectral avar calculation
#   
#   ##calculate tapers for this data spacing
#   V.mat <- get_tapers(t.vec, W = myW, K = numTapers) #W was 4/N
#   
#   MTSE_full <- MT_spectralEstimate_fft(x.t, freq, V.mat$tapers) 
#   
#   Cov.mat=NA
#   ### calculate the covariance matrix 
#   if(calcCov==T){
#     Cov.mat=spec_cov.mat(x.t, t.vec, N.fourier, V.mat$tapers, isWhite,acf.lag,numCores)
#     # Cov.mat=spec_cov.mat_slow(x.t, t.vec, N.fourier, V.mat$tapers, isWhite,acf.lag)
#   }
#   
#   return(list(freq=freq,
#               spec.hat=MTSE_full$spectrum,Cov.mat=Cov.mat))
#   
# }
# 
# 
# spec_cov.mat_slow <- function(X.t, t.vec, N.fourier, taperMat, isWhite,acf.lag){
#   N <- length(stats::na.omit(X.t)) #length of data without gaps
#   K <- dim(taperMat)[2] 
#   
#   freq <- seq(0, 0.5, length.out = N.fourier)
#   
#   C.mat <- matrix(NA, nrow = N.fourier, ncol = N.fourier)
#   
#   if(isWhite){
#     R_mat <- diag(1, nrow = N) #to start
#   }
#   if(!isWhite){
#     s_acf <- stats::acf(X.t, plot=FALSE, lag.max=acf.lag,na.action = stats::na.exclude)$acf
#     
#     # Create a Toeplitz matrix from the autocorrelation values
#     R_mat <- matrix(0, nrow = N, ncol = N)
#     R_mat <- stats::toeplitz(c(s_acf, rep(0, N - acf.lag - 1)))
#   }
#   
#   for(i in 1:N.fourier){
#     if(i %% 100 == 0){print(paste(i," of ",N.fourier))}
#     j = 1
#     while(j <= i){
#       V_star_mat <- Conj(t(taperMat*exp(-1i*2*pi*freq[i]*t.vec)))
#       V_mat <- taperMat*exp(-1i*2*pi*freq[j]*t.vec)
#       C.mat[i,j] <- (1/K)*norm(V_star_mat%*%R_mat%*%V_mat, type = "2") 
#       j = j+1
#     }
#   }
#   C.mat[upper.tri(C.mat)] <- t(C.mat)[upper.tri(C.mat)]
#   
#   return(C.mat = C.mat)
# }
# 
# 
# 
# #########################################################
# ### start of new parallel cov Functions
# #########################################################
# 
# ## upper triangle indices
# generate_upper_triangle_indices <- function(N) {
#   
#   print("generate_upper_triangle_indices")
#   # Create a sequence of row indices
#   row_indices <- rep(1:N, each = N)
#   # Create a sequence of column indices
#   col_indices <- rep(1:N, times = N)
#   # Combine row and column indices into a matrix
#   all_indices <- cbind(row_indices, col_indices)
#   # Filter out indices to get only the upper triangle (excluding the diagonal)
#   upper_triangle_indices <- all_indices[row_indices < col_indices, ]
#   return(upper_triangle_indices)
# }
# 
# ## compute the entry C.mat_ij
# compute_entry_parallel <- function(ij, taperMatrix, sett.vec=t.vec, setfreq=freq,setK = K, setN = N, c = c_vec){
#   #i,j, taperMat = taperMatrix, setK = K, setN = N, c = c_vec
#   print("compute_entry_parallel")
#   
#   i = ij[1]
#   j = ij[2]
#   
#   V_star <- t(taperMatrix*exp(1i*2*pi*setfreq[i]*sett.vec))
#   V <- taperMatrix*exp(-1i*2*pi*setfreq[j]*sett.vec)
#   return(est_entry_FFT(V_star, V, c, setK, setN))
# }
# 
# ## fill in the upper triangle of a matrix from a vector
# fill_upper_triangle <- function(vec, N, incl_diag = TRUE) {
#   
#   print("fill_upper_triangle")
#   # Initialize an NxN matrix with NA
#   mat <- matrix(NA, nrow = N, ncol = N)
#   
#   # Get the indices of the upper triangle (excluding the diagonal)
#   upper_tri_indices1 <- which(upper.tri(mat, diag = incl_diag), arr.ind = TRUE)
#   upper_tri_indices <- upper_tri_indices1[order(upper_tri_indices1[,1]),]
#   
#   # Check if the length of the vector matches the number of upper triangle elements
#   if (length(vec) != nrow(upper_tri_indices)) {
#     stop("Length of vector does not match the number of upper triangle elements excluding the diagonal.")
#   }
#   
#   # Fill the matrix upper triangle with the elements of the vector
#   mat[upper_tri_indices] <- vec
#   
#   return(mat)
# }
# 
# spec_cov.mat=function(x.t, t.vec,N.fourier,taperMat,isWhite,
#                       max.lag.acf,numCores){
#   
# 
#   print("spec_cov.mat")
#   N <- length(stats::na.omit(x.t)) #length of data without gaps
#   taperMatrix=taperMat
#   K <- dim(taperMatrix)[2] 
#   
#   freq <- seq(0, 0.5, length.out = N.fourier)
#   
#   # t.vec <- 1:length(x.t)
#   # t.vec[which(is.na(x.t))] <- NA
#   # t.vec <- stats::na.omit(t.vec)
#   # numCores <- numCores
#   
#   ### make the cluster
#   cl <- parallel::makeCluster(numCores)
#   
#   ## Pre-processing ####
#   
#   ### 1. load packages ####
#   parallel::clusterEvalQ(cl, {
#     library(Rcpp)
#   })
#   
#   ### 2. send source C code ####
#   parallel::clusterEvalQ(cl, {
#     Rcpp::sourceCpp('/home/aak3/NIST/atomic-clock/CovarianceCalculation/est_entry_FFT.cpp')
#     # Rcpp::sourceCpp(paste(folderLocation,'CovarianceCalculation/est_entry_FFT.cpp',sep=""))
#   })
#   
#   ### 3. Calculate predetermined variables ####
#   
#   #### c vector
#   print("c vector")
#   if(isWhite){
#     sample.acf <- c(1,rep(0,max.lag.acf-1))
#   }
#   if(!isWhite){
#     sample.acf <- stats::acf(x.t, plot=FALSE, lag.max=max.lag.acf,na.action = stats::na.exclude)$acf
#   }
#   
#   c_vec <- c(sample.acf,rep(0, times = N-max.lag.acf), 0, rep(0, times = N-max.lag.acf), rev(sample.acf[-1]))
#   
#   #### list of indices (don't really need to send this to the cluster, but good for debugging if needed)
#   upper_triangle_indices <- generate_upper_triangle_indices(N.fourier)
#   
#   ### 4. Send variables/necessary functions to the cluster ####
#   print("K is")
#   print(K)
#   print("N is")
#   print(N)
#   parallel::clusterExport(cl, list("t.vec", "N", "K", "c_vec", "taperMatrix", "freq", 
#                                    "upper_triangle_indices"),envir=environment()) #objects
#   
#   ######################this is not exporting correctly!! adding taperMatrix explicitly below seemed to help, but then freq was the holdup
#   # i feel like maybe the variables aren't talking to the outside functions, since i make it all the way to "  print("start parallel")" without error
#   # i tried adding the functions inside this function, did not help
#   # maybe try pass all variables as arguments of your function CrossValJM() so that they are automatically exported to the clusters
#   # https://stackoverflow.com/questions/47424367/error-in-check-for-remote-errors-val-5-nodes-produced-an-error-object-not-fo
#   print("exported the objects")
#   parallel::clusterExport(cl, "compute_entry_parallel") #functions
# 
#   ## Run code on cluster ####
#   
#   ##### (run this code chunk all at once ###
#   #####   for most accurate timing)      ###
#   start_time_fast = Sys.time() #for timing, start time
#   print("start parallel")
#   my_list <- list() #a vector to hold the entries
#   my_list <- parallel::parApply(cl, FUN = compute_entry_parallel, X = upper_triangle_indices, MARGIN = 1)#,taperMatrix=taperMat) this helped, but then got stuck at freq, not a good long term fix
#   
#   total_time_fast = Sys.time() - start_time_fast
#   print(total_time_fast) 
#   
#   parallel::stopCluster(cl) #stop the cluster
#   #######--------------------------------###
#   
#   
#   
#   ## Make the Covariance matrix ####
#   ### 1. Fill in upper triangle of C.mat ####
#   C.mat <- fill_upper_triangle(vec = my_list, N.fourier, incl_diag = FALSE)
#   
#   ### 2. Fill in lower triangle of C.mat ####
#   C.mat[lower.tri(C.mat)] <- t(C.mat)[lower.tri(C.mat)]
#   
#   ### 3. Fill in Diagonal of C.mat ####
#   tN = max.lag.acf
#   R_mat <- stats::toeplitz(c(seq(1,0.1, length.out = tN), rep(0, times = N-tN))) #to start
#   diag(C.mat) <- norm(t(taperMatrix)%*%R_mat%*%taperMatrix/K, type = "2")
#   
#   ## look at result ##
#   return(list(C.mat = C.mat))
# }
# #########################################################
# ### end of new parallel cov Functions
# #########################################################
# 

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


spec_cov.mat_slow_WN <- function(t.vec, N.fourier, taperMat){
  K <- dim(taperMat)[2]
  
  freq <- seq(0, 0.5, length.out = N.fourier)
  
  C.mat <- matrix(NA, nrow = N.fourier, ncol = N.fourier)
  
  for(i in 1:N.fourier){
    if(i %% 100 == 0){print(paste(i," of ",N.fourier))}
    j = 1
    while(j <= i){
      V_star_mat <- Conj(t(taperMat*exp(-1i*2*pi*freq[i]*t.vec)))
      V_mat <- taperMat*exp(-1i*2*pi*freq[j]*t.vec)
      
      V_star_matOMIT=V_star_mat[,-which(colSums(is.na(V_star_mat))==5)]
      V_matOMIT=V_mat[-which(rowSums(is.na(V_mat))==5),]
      
      C.mat[i,j] <- (1/K)*norm(V_star_matOMIT%*%V_matOMIT, type = "2") 
      j = j+1
    }
  }
  C.mat[upper.tri(C.mat)] <- t(C.mat)[upper.tri(C.mat)]
  
  return(C.mat = C.mat)
}



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
    CI.limits <- dplyr::bind_rows("lower" = s.2 - s.2/N, "upper" = s.2 +  s.2/N)
  }else if(avar_type == "chisquared"){
    CI.limits <- dplyr::bind_rows("lower" = s.2*edf/stats::qchisq(1-a,edf),"upper" = s.2*edf/stats::qchisq(a, edf) )
  }else{
    warning("Invalid avar_type, should be 'simple' or 'chisquared'")
  }
  
  CI.limits = dplyr::bind_cols(data.frame(tau=taus),CI.limits)
  
  return(CI.limits)
}
