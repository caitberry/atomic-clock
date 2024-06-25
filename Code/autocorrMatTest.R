
##############corr mat

# Generate example time series data
set.seed(123)
ts_data <- ts(rnorm(10000))  # Replace with your time series data

# Calculate autocorrelations up to lag p
p <- 4  # Define the maximum lag
acf_result <- acf(ts_data, lag.max = p, plot = FALSE)

# Extract autocorrelation values
acf_values <- acf_result$acf

# Get the length of the time series data
N <- length(ts_data)

# # Initialize the autocorrelation matrix with zeros
autocorr_matrix <- matrix(0, nrow = N, ncol = N)

# Fill in the autocorrelation values for the specified lags
# 
# for (lag in 0:p) {
#   for (i in 1:(N - lag)) {
#     autocorr_matrix[i, i + lag] <- acf_values[lag + 1]
#     autocorr_matrix[i + lag, i] <- acf_values[lag + 1]  # Since the matrix is symmetric
#   }
# }
# 
# autocorr_matrix-toeplitz_matrix

# Create a Toeplitz matrix from the autocorrelation values
autocorr_matrix <- toeplitz(c(acf_values, rep(0, N - p - 1)))

head(autocorr_matrix)
diag(autocorr_matrix)
autocorr_matrix[1:10,1:10]


####### try with real data

###2024-04-26-ErYb-AlSr

dataFile="2024-04-26-ErYb-AlSr"
allDat=read.csv(paste("/home/aak3/NIST/atomic-clock/Data/clockRatioOffsetTimeseries_first2024dryRun/",dataFile,".csv",sep=""),header = F,sep="",na.strings = "NaN",skip = 3)
# need to create dat vector with NAs to send to 2.

### 03/05/2018 and 03/06/2018 analyzed on 07/10/18 frac.diff. = (measured-ideal)/ideal * 1e15
dat=allDat$V2

# Generate example time series data
ts_data <- ts(dat)  # Replace with your time series data

# Calculate autocorrelations up to lag p
p <- 4  # Define the maximum lag
acf_result <- acf(ts_data, lag.max = p, plot = FALSE)

# Extract autocorrelation values
acf_values <- acf_result$acf

# Get the length of the time series data
N <- length(ts_data)

# # Initialize the autocorrelation matrix with zeros
autocorr_matrix <- matrix(0, nrow = N, ncol = N)

# Create a Toeplitz matrix from the autocorrelation values
autocorr_matrix <- toeplitz(c(acf_values, rep(0, N - p - 1)))

head(autocorr_matrix)
diag(autocorr_matrix)
autocorr_matrix[1:10,1:10]