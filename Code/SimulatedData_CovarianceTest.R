set.seed(123)
N=1000 #length of data with gaps

x.t <- rnorm(N)  # Replace with your time series data
omitted=c(301:400)
x.t[omitted]=NA

t.vec <- 1:N
t.vec[which(is.na(x.t))] <- NA

mtse=modules::use("/home/aak3/NIST/atomic-clock/Functions.R")
source("/home/aak3/NIST/atomic-clock/Functions_SpectrumCovariance.R")

W=4/N*3
K=5
V.mat <- mtse$get_tapers(t.vec, W = W, K = K) 
print(V.mat$e.values)
taperMatrix=V.mat$tapers

MTSE_full <- mtse$MT_spectralEstimate_fft(x.t, taperMatrix) 
plot(MTSE_full$freqs,MTSE_full$spectrum)

Cov.mat=spec_cov.mat(x.t = x.t, t.vec = t.vec, N.fourier = floor(N/2) + 1, 
                     taperMat = V.mat$tapers, isWhite = isWhite,
                     max.lag.acf = 4,
                     numCores=80)

