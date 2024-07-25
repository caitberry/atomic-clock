############################## Parallelization ##################################
##### Title: Parallel_Covariance_Windows.R                                  #####
##### Description: This file contains code to parallelize the calculation   #####
#####              of the covariance matrix of the spectral estimate using  #####
#####              c++ function and parApply() in the "parallel" package.   #####
#####              This is meant to be used on a windows machine as an      #####
#####              alternative to mcapply() which doesn't work on Windows.  #####
#################################################################################-
#Set working directory for Cait: 
# setwd("/home/cmb15/ClockDataAnalysis/Code/Cait")
# source("Parallel_Covariance_Windows.R")

###################-
#### libraries ####
###################-

library(parallel)
library(Rcpp)
library(RcppArmadillo)

###################-
#### functions ####
###################-

## upper triangle indices
generate_upper_triangle_indices <- function(N) {
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
  
  V_star <- t(taperMat*exp(1i*2*pi*freq[i]*t.vec))
  V <- taperMat*exp(-1i*2*pi*freq[j]*t.vec)
  return(est_entry_FFT(V_star, V, c, setK, setN))
}

## fill in the upper triangle of a matrix from a vector
fill_upper_triangle <- function(vec, N, incl_diag = TRUE) {
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

## mtse module
mtse <- modules::use("Functions.R")

###################-
#####  Method #####
###################-


## Start a Cluster ####

### pick number of cores 
### can use detectCores() to see how many are available
numCores <- 14

### make the cluster
cl <- makeCluster(numCores)

## Pre-processing ####

### 1. load packages ####
clusterEvalQ(cl, {
  library(Rcpp)
})

### 2. send source C code ####
clusterEvalQ(cl, {
  sourceCpp('est_entry_FFT.cpp')
})

### 3. Calculate predetermined variables ####

N=128 # length of data
N.fourier = floor(N/2) + 1
t.vec <- 1:N  # time vector

#### tapers
setW = 8/N
K = 5
taperObject <- mtse$get_tapers(t.vec, W = setW, K = K)
taperMatrix <- taperObject$tapers

#### c vector
max.lag.acf = 10
sample.acf <- seq(1,0.1, length.out = max.lag.acf) #toy example
c_vec <- c(sample.acf,rep(0, times = N-max.lag.acf), 0, rep(0, times = N-max.lag.acf), rev(sample.acf[-1]))
length(c_vec)

#### frequency vector
freq <- seq(0, 0.5, length.out = N.fourier)  # Frequency vector

#### list of indices (don't really need to send this to the cluster, but good for debugging if needed)
upper_triangle_indices <- generate_upper_triangle_indices(N.fourier)

### 4. Send variables/necessary functions to the cluster ####

clusterExport(cl, list("t.vec", "N", "K", "c_vec", "taperMatrix", "freq", "upper_triangle_indices")) #objects
clusterExport(cl, "compute_entry_parallel") #functions

start_time_fast = Sys.time()
compute_entry_parallel(upper_triangle_indices[1,])
total_time_fast = Sys.time() - start_time_fast
total_time_fast


## Run code on cluster ####

##### (run this code chunk all at once ###
#####   for most accurate timing)      ###
  start_time_fast = Sys.time() #for timing, start time
  print("start parallel")
  my_list <- list() #a vector to hold the entries
  my_list <- parApply(cl, FUN = compute_entry_parallel, X = upper_triangle_indices, MARGIN = 1)
  
  total_time_fast = Sys.time() - start_time_fast
  print(total_time_fast) 
  
  stopCluster(cl) #stop the cluster
#######--------------------------------###
  

  
## Make the Covariance matrix ####
### 1. Fill in upper triangle of C.mat ####
C.mat <- fill_upper_triangle(vec = my_list, N.fourier, incl_diag = FALSE)

### 2. Fill in lower triangle of C.mat ####
C.mat[lower.tri(C.mat)] <- t(C.mat)[lower.tri(C.mat)]

### 3. Fill in Diagonal of C.mat ####
tN = max.lag.acf
R_mat <- toeplitz(c(seq(1,0.1, length.out = tN), rep(0, times = N-tN))) #to start
diag(C.mat) <- norm(t(taperMatrix)%*%R_mat%*%taperMatrix/K, type = "2")

## look at result ##
C.mat[1:10,1:10]
