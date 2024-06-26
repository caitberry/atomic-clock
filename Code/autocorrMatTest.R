
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

dat=smushSmallGaps(dat)
# Generate example time series data
ts_data <- c(dat,rep(NA,10))  # Replace with your time series data

# Calculate autocorrelations up to lag p
p <- 4  # Define the maximum lag
acf_result <- acf(ts_data, lag.max = p, plot = FALSE,na.action = na.pass)

# find gaps

# Function to split a time series into chunks based on NAs
split_ts_by_na <- function(ts_data) {
  # Find the positions of NAs
  na_positions <- which(is.na(ts_data))

  na_positions = na_positions[which(diff(na_positions)>1)]
  # Initialize start and end points of chunks
  chunk_starts <- c(1, na_positions + 1)
  chunk_ends <- c(na_positions - 1, length(ts_data))
  
  # Remove invalid start and end positions
  chunk_starts <- chunk_starts[chunk_starts <= length(ts_data)]
  chunk_ends <- chunk_ends[chunk_ends >= 1]
  
  # Ensure each start position is less than or equal to the corresponding end position
  valid_chunks <- chunk_starts <= chunk_ends
  chunk_starts <- chunk_starts[valid_chunks]
  chunk_ends <- chunk_ends[valid_chunks]
  
  # Split the time series into chunks
  chunks <- mapply(function(start, end) ts_data[start:end], chunk_starts, chunk_ends, SIMPLIFY = FALSE)
  
  # Remove empty chunks
  chunks <- Filter(function(x) length(x) > 0, chunks)
  
  return(chunks)
}

# # Generate example time series data with NAs
# set.seed(123)
# ts_data <- ts(c(NA, NA, 1:20, rep(NA,10), 21:30, NA, NA, 31:44, rep(NA,10), 45:77, NA, NA, NA))

# Apply the function to the example data
chunks <- split_ts_by_na(ts_data)

# Print the chunks
print(chunks)


diff(which(is.na(chunks[[1]])))
diff(which(is.na(chunks[[2]])))
diff(which(is.na(chunks[[3]])))
diff(which(is.na(chunks[[4]])))
diff(which(is.na(chunks[[5]])))

plot(ts_data)
plot(chunks[[2]])
plot(chunks[[3]])
plot(chunks[[4]])
plot(chunks[[5]])
tail(ts_data)

# Calculate autocorrelations up to lag p
p <- 4  # Define the maximum lag

acf1=acf(chunks[[2]],na.action = na.omit,lag.max = p,plot = F)
acf2=acf(chunks[[3]],na.action = na.omit,lag.max = p,plot = F)
acf3=acf(chunks[[4]],na.action = na.omit,lag.max = p,plot = F)
acf4=acf(chunks[[5]],na.action = na.omit,lag.max = p,plot = F)

acf1$lag
as.numeric(test)
acfOut=data.frame(chunk=rep(1:4,each=p+1),
                  acf=c(as.numeric(acf1$acf),
                        as.numeric(acf2$acf),
                        as.numeric(acf3$acf),
                        as.numeric(acf4$acf)),
                  lag=c(as.numeric(acf1$lag)),
                  N=)

ggplot(acfOut,aes(lag,acf,col=chunk))+
  geom_point()

test=acf(chunks[[2]],na.action = na.omit,lag.max = p,plot = F)
test$n.used






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
