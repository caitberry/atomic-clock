
library(parallel)
library(Rcpp)
library(RcppArmadillo)

## Spectral Estimate and Uncertainty ################################################
#inputs:  (->) x.t = vector of data (possibly with NA values, include if they are in the dataset)
#         (->) t.vec = vector of data (possibly with NA values, include if they are in the dataset)
#         (->) N.fourier
#         (->) numTapers = number of tapers
#         (->) calcCov = True if you want to calculate the covariance matrix
#         (->) myW = Analysis half-bandwidth
#         (->) isWhite
#         (->) acf.lag
#output:  (<-) spectral estimate with covariance matrix

spectralEstWithUnc <- function(x.t,t.vec,N.fourier,numTapers,calcCov=T,myW,isWhite = TRUE,acf.lag=4,numCores){
  
  freq <- seq(0,0.5, length.out = N.fourier)
  
  ##calculate tapers for this data spacing
  V.mat <- mtse$get_tapers(t.vec, W = myW, K = numTapers) #W was 4/N
  
  MTSE_full <- mtse$MT_spectralEstimate_fft(x.t, V.mat$tapers) 
  
  Cov.mat=NA
  ### calculate the covariance matrix 
  if(calcCov==T){
    Cov.mat=spec_cov.mat(x.t, t.vec, N.fourier, V.mat$tapers, isWhite,acf.lag,numCores)
  }
  
  return(list(freq=freq,
              spec.hat=MTSE_full$spectrum,Cov.mat=Cov.mat, e.values=V.mat$e.values))
}

#########################################################
### start of new parallel cov Functions
#########################################################

## upper triangle indices
generate_upper_triangle_indices <- function(N) {
  
  # print("generate_upper_triangle_indices")
  # Create a sequence of row indices
  row_indices <- rep(1:N, each = N)
  # Create a sequence of column indices
  col_indices <- rep(1:N, times = N)
  # Combine row and column indices into a matrix
  all_indices <- cbind(row_indices, col_indices)
  # Filter out indices to get only the upper triangle (excluding the diagonal)
  upper_triangle_indices <- all_indices[row_indices < col_indices, ]
  return(upper_triangle_indices)
}

## compute the entry C.mat_ij
compute_entry_parallel <- function(ij, taperMat = taperMatrix, setK = K, setN = N, c = c_vec){
  #i,j, taperMat = taperMatrix, setK = K, setN = N, c = c_vec
  i = ij[1]
  j = ij[2]
  
  N.noNA=sum(!is.na(t.vec))
  taperMat=na.omit(taperMat)
  t.vec=na.omit(t.vec)
  
  U_star <- t(taperMat*exp(1i*2*pi*freq[i]*t.vec))
  U <- taperMat*exp(-1i*2*pi*freq[i]*t.vec)
  V_star <- t(taperMat*exp(1i*2*pi*freq[j]*t.vec))
  V <- taperMat*exp(-1i*2*pi*freq[j]*t.vec)
  
  ### NEW CALC needs to go in est_entry_FFT
  
  return(est_entry_FFT(V_star, V, c, setK, N.noNA))
}


## fill in the upper triangle of a matrix from a vector
fill_upper_triangle <- function(vec, N, incl_diag = TRUE) {
  
  # print("fill_upper_triangle")
  # Initialize an NxN matrix with NA
  mat <- matrix(NA, nrow = N, ncol = N)
  
  # Get the indices of the upper triangle (excluding the diagonal)
  upper_tri_indices1 <- which(upper.tri(mat, diag = incl_diag), arr.ind = TRUE)
  upper_tri_indices <- upper_tri_indices1[order(upper_tri_indices1[,1]),]
  
  # Check if the length of the vector matches the number of upper triangle elements
  if (length(vec) != nrow(upper_tri_indices)) {
    stop("Length of vector does not match the number of upper triangle elements excluding the diagonal.")
  }
  
  # Fill the matrix upper triangle with the elements of the vector
  mat[upper_tri_indices] <- vec
  
  return(mat)
}

spec_cov.mat=function(x.t, t.vec,N.fourier,taperMat,isWhite,
                      max.lag.acf,numCores){
  
  
  # print("spec_cov.mat")
  N <- length(x.t) #length of data with gaps
  taperMatrix=taperMat
  K <- dim(taperMatrix)[2] 
  
  freq <- seq(0, 0.5, length.out = N.fourier)
  
  ### make the cluster
  cl <- parallel::makeCluster(numCores)
  
  ## Pre-processing ####
  
  ### 1. load packages ####
  parallel::clusterEvalQ(cl, {
    library(Rcpp)
  })
  
  ### 2. send source C code ####
  parallel::clusterEvalQ(cl, {
    Rcpp::sourceCpp('CovarianceCalculation/est_entry_FFT.cpp')
    # Rcpp::sourceCpp(paste(folderLocation,'CovarianceCalculation/est_entry_FFT.cpp',sep=""))
  })
  
  ### 3. Calculate predetermined variables ####
  
  #### c vector
  # print("c vector")
  if(isWhite){
    sample.acf <- c(1,rep(0,max.lag.acf-1))
  }
  if(!isWhite){
    sample.acf <- stats::acf(x.t, plot=FALSE, lag.max=max.lag.acf,na.action = stats::na.exclude)$acf
  }
  
  N.noNA=sum(!is.na(t.vec))
  c_vec <- c(sample.acf,rep(0, times = N.noNA-max.lag.acf-1), 0, rep(0, times = N.noNA-max.lag.acf-1), rev(sample.acf[-1]))
  
  #### list of indices (don't really need to send this to the cluster, but good for debugging if needed)
  upper_triangle_indices <- generate_upper_triangle_indices(N.fourier)
  
  ### 4. Send variables/necessary functions to the cluster ####
  print("K is")
  print(K)
  print("N is")
  print(N)
  parallel::clusterExport(cl, list("t.vec", "N", "K", "c_vec", "taperMatrix", "freq", 
                                   "upper_triangle_indices"),envir=environment()) #objects
  
  print("exported the objects")
  parallel::clusterExport(cl, "compute_entry_parallel") #functions
  
  ## Run code on cluster ####
  
  ##### (run this code chunk all at once ###
  #####   for most accurate timing)      ###
  start_time_fast = Sys.time() #for timing, start time
  print("start parallel")
  my_list <- list() #a vector to hold the entries
  my_list <- parallel::parApply(cl, FUN = compute_entry_parallel, X = upper_triangle_indices, MARGIN = 1)
  
  total_time_fast = Sys.time() - start_time_fast
  print(total_time_fast) 
  
  parallel::stopCluster(cl) #stop the cluster
  #######--------------------------------###
  
  
  
  ## Make the Covariance matrix ####
  ### 1. Fill in upper triangle of C.mat ####
  C.mat <- fill_upper_triangle(vec = my_list, N.fourier, incl_diag = FALSE)
  
  ### 2. Fill in lower triangle of C.mat ####
  C.mat[lower.tri(C.mat)] <- t(C.mat)[lower.tri(C.mat)]
  
  ### 3. Fill in Diagonal of C.mat ####
  # tN = max.lag.acf
  # R_mat <- toeplitz(c(seq(1,0.1, length.out = tN), rep(0, times = N-tN))) #to start
  s_acf <- stats::acf(x.t, plot=FALSE, lag.max=max.lag.acf,na.action = stats::na.exclude)$acf
  
  # Create a Toeplitz matrix from the autocorrelation values
  R_mat <- matrix(0, nrow = N.noNA, ncol = N.noNA)
  R_mat <- stats::toeplitz(c(s_acf, rep(0, N.noNA - max.lag.acf - 1)))

  # Calculate variances
  taperMat=stats::na.omit(taperMatrix)
  t.vec=stats::na.omit(t.vec)
  
  for (j in 1:length(j)){
    V_star <- t(taperMat*exp(1i*2*pi*freq[j]*t.vec))
    V <- taperMat*exp(-1i*2*pi*freq[j]*t.vec)
    
    diag(C.mat)[j] <- (1/K)*2*tr(R_mat %*% V %*% V_star %*%  V %*% V_star %*% R_mat)
  } #check if this is different for diff frequencies, if not just calc one
  
  ## look at result ##
  return(C.mat)
}
#########################################################
### end of new parallel cov Functions
#########################################################

