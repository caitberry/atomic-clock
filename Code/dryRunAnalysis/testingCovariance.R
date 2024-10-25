
folderLocation="C:/Users/aak3/Documents/atomic-clock/"
# folderLocation="/home/aak3/NIST/atomic-clock/"

# source(file = paste(folderLocation,"Functions_SpectrumCovariance.R",sep="")) 

source(file = paste(folderLocation,"Code/analysisFunctions.R",sep=""))

### data analysis steps
# 1. read in data

N=500
x.t=rnorm(N)
t.vec=1:N

plot(t.vec,x.t)

########################################################################################
# 2. concatenate small gaps
########################################################################################

# formattedData=formatDat(dat)
# 
# x.t=formattedData$x.t
# t.vec=formattedData$t.vec
# N=formattedData$N
# 
# plot(t.vec,x.t)
# 
########################################################################################
# 3. get spectral estimate
# 3a) try a variety of W and K values and check eigenvalues to make sure you've made good selections
# can start with 8/N or 12/N for W
########################################################################################
myW=4/N*3
myK=10
V.mat <- mtse$get_tapers(t.vec, W = myW, K = myK) 
V.mat$e.values
startTime=Sys.time()

# myDatSpecEst=spectralEstWithUnc(x.t = x.t,t.vec=t.vec,N.fourier = floor(N/2) + 1,
#                        numTapers = myK,calcCov = T,
#                        myW = myW,isWhite = F,acf.lag = 4,numCores=80)
# 
#############################covariance calc

#### tapers
setW = myW
K = myK
taperMatrix=V.mat$tapers

#### frequency vector
N.fourier = floor(N/2) + 1
freq <- seq(0, 0.5, length.out = N.fourier)  # Frequency vector

### pick number of cores
### can use detectCores() to see how many are available
numCores <- 10

# isWhite=F
max.lag.acf=4
sample.acf <- stats::acf(x.t, plot=FALSE, lag.max=max.lag.acf,na.action = stats::na.exclude)$acf

source(paste(folderLocation,"CovarianceCalculation/Parallel_Covariance_Windows.R",sep=""))

Sys.time()-startTime
