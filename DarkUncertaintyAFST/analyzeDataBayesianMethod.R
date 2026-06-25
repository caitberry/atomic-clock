############################################################################
### 
### Dark uncertainty for atomic clock data
### Analyze simulated data using Bayesian method
### 
### Angela Folz
### April 2026
### Code adapted from Amanda Koepke
### 
############################################################################

library(readr)
library(dplyr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

path <- "/Users/adf2/OneDrive - NIST/Documents/SpectralAnalysis/atomic-clock/"
gitpath <- "Code/ComparisonAnalysis2025/DarkUncertaintyMethods/"
simdatfolder <- "DarkUncertaintyAFST/simulatedData/"
successmetricsfolder <- "DarkUncertaintyAFST/successMetrics/"
BayesCIfolder <- "DarkUncertaintyAFST/CIouts/"

todayDate <- format(Sys.Date(), "%Y%m%d")

############################################################################
### Load simulated data from random effects model
############################################################################

# choose data to analyze (which model it was simulated under) (ONE)
howDatSim <- "RandomEffects"
# howDatSim <- "MulBirge"

simfiles <- list.files(path=paste0(path, simdatfolder), pattern=c(howDatSim,".csv")) # get file paths
simdat <- lapply(paste0(path, simdatfolder, simfiles), ifelse(howDatSim=="MulBirge", read.csv2, read.csv))
names(simdat) <- basename(simfiles) # name the list elements after the files
# simdat <- read_csv(paste0(path, simdatfolder, simfiles), id="source_file") # read and combine into one large df

n_iter <- 10 # 100 #10000 # simdat[[1]][nrow(simdat[[1]]),"iter"] # get number of simulation iterations #!!!!!!!!!!!!!!!!!!!

############################################################################
### Parameters
############################################################################

analysisMethod <- "Bayes"

sampling.iter <- 10000 # 5000 # NOTE: Amanda's numberOfSims = my n_iter; Amanda's num.iters is for iter in fitGaussGauss

# every combo of N and xi
allparams <- expand.grid(N_set = unique(sapply(simdat, function(x) x[["N"]][1])), # get each N from simdat
                         xi_set = unique(sapply(simdat, function(x) x[["xi"]][1]))) # get each xi from simdat

myparams <- allparams %>%
  mutate(adaptDeltaVals = #########right now, setting these based on N # TODO: good values??????????????????????????
           ifelse(N_set == 13, .95, 
                  ifelse(N_set== 33, .8,
                         .8))) %>%
  mutate(maxTreeDepthVals = #########right now, setting these based on N
           ifelse(N_set== 13, 15,
                  ifelse(N_set== 33, 10,
                         10))) %>%
  mutate(xidf_set = 
           ifelse(N_set == 5 & xi_set %in% c(1,2), 2,
                  4)) #set prior sd to 4 for all except small N and small xi case, then too big
myparams <- myparams[c(1,4,7,2,5,8,3,6,9),] # reorder to match sidmdat


# TODO: check if need these !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ratiolab = "YbSr"
lblabel_set = paste(ratiolab,"_minunc",sep="")
ublabel_set = paste(ratiolab,"_maxunc",sep="")
alpha_set = NA
nu_set = NA

############################################################################
### Functions
############################################################################

# fitGaussGauss
source(file = paste(path,"DarkUncertaintyAFST/fitBayesianModels.R", sep=""))

#################################################################

# wrapper to catch errors and keep going
oneTrial <- function(dat, GGmodel, N, xi, mu, num.iters, alpha, nu, 
                     howDatSim, adaptDeltaVal, maxTreeDepthVal, xidf, idx){
  tmp_df <- tryCatch(runOneTrial(dat, GGmodel, N, xi, mu, num.iters, alpha, nu, 
                                 howDatSim, adaptDeltaVal, maxTreeDepthVal, xidf, idx), 
                     error = function(err){print("err"); return(data.frame("method"=NA))})
  return(tmp_df)
}

#################################################################

runOneTrial <- function(dat, GGmodel, N, xi, mu, num.iters, alpha, nu, 
                        howDatSim, adaptDeltaVal, maxTreeDepthVal, xidf, idx){
  
  stan_data <- list(N = N,
                    x = dat$x,
                    u = dat$u,
                    tdf = xidf) # note, xi is called tau in fitGaussGauss
  # print(stan_data)
  GGres <- fitGaussGauss(stanDat = stan_data, trueTau = xi, trueMu = mu, model = GGmodel, n.iter = num.iters,
                         adapt.delta = adaptDeltaVal, tree.depth = maxTreeDepthVal, N.obs = N, iter.idx = idx)
  bayesRes <- GGres
  bayesRes[bayesRes$parameter == "tau", "parameter"] <- "xi" # rename tau to xi
  
  newOut <- bayesRes
  newOut$simulation <- howDatSim
  newOut$N <- N
  newOut$mu <- mu
  newOut$xi <- xi
  newOut$num.iters <- num.iters
  newOut$alpha <- alpha
  newOut$nu <- nu
  newOut$xidf <- xidf
  newOut$maxtree <- maxTreeDepthVal
  newOut$adaptDelta <- adaptDeltaVal
  
  return(newOut)
}

#################################################################

runForManyRepeats <- function(dat, numberOfSims, num.iters, alpha, nu, howDatSim, 
                              theadaptdelta, thetreedepth, thexidf, today, lblabel, ublabel){
  startTime=Sys.time()
  
  N = dat[[1]]$N[1]
  xi = dat[[1]]$xi[1]
  mu = dat[[1]]$mu[1]
  lb = dat[[1]]$lb[1]
  ub = dat[[1]]$ub[1]
  
  print(paste("xi =",xi, ", N =",N,", adapt delta =",theadaptdelta,
              ", max tree =",thetreedepth,", ub =",ub,", lb =",lb))
  
  myHGG_noDoF_reparam <- stan_model(paste0(path, "DarkUncertaintyAFST/", "myHGG_noDoF_reparam.stan"))
  
  CIout <- list()
  
  for(i in 1:numberOfSims){

    CIout[[i]] <- oneTrial(dat = dat[[1]][((i-1)*N+1):(i*N), c("Day", "x", "u")],
                           GGmodel = myHGG_noDoF_reparam,
                           N = N, xi = xi, mu = mu, num.iters = num.iters, alpha = alpha, nu = nu,
                           howDatSim = howDatSim,
                           adaptDeltaVal = theadaptdelta, maxTreeDepthVal = thetreedepth, xidf = thexidf, idx=i)

    CIout[[i]]$Iteration <- i
    CIout[[i]]$lb <- lb
    CIout[[i]]$ub <- ub

    if(i %% 10==0) { # save 10 at a time
      print(paste("iteration", i))
      saveCIout <- bind_rows(CIout[(i-9):i])
      write_csv(saveCIout, paste(path, BayesCIfolder, howDatSim, "CIout", today, "_xi",xi, "_N",N, "_",n_iter,"iter",
                                "_ub",ublabel,"_lb",lblabel, ".csv", sep=""), append=T,
                col_names=!file.exists(paste(path, BayesCIfolder, howDatSim, "CIout", today, "_xi",xi, "_N",N,
                                             "_",n_iter,"iter", "_ub",ublabel,"_lb",lblabel, ".csv", sep="")))
    }
  }
  
  runTime <- round(difftime(Sys.time(), startTime, units = "mins")[[1]], 3)
  print(paste("run time =", runTime, "min"))
  cat(paste("N = ",N, ", xi = ",xi, ", ",n_iter," iterations: run time = ", runTime, " min\n", sep=""),
      file = paste(path, BayesCIfolder, howDatSim, "times", today, ".txt", sep=""), append=T)
  
  return("Done!")
  
}


############################################################################
### Analyze simulated data using Bayesian method
############################################################################

##########################################
# Test whether method works

# # one iteration, 1 factor combo
# 
# # choose data set (parameter combination) to work with
# dat <- as.data.frame(simdat[simfiles[6]])
# names(dat) <- c("Day", "x", "u", "N", "mu", "xi", "lb","ub","seed", "iter") # fix names
# 
# it <- 1 # choose which iteration
# dat_1iter <- dat[dat$iter==it,]
# 
# outN13xi3_1iter_bayes <- runForManyRepeats(dat_1iter,
#                                      numberOfSims = 1,
#                                      num.iters = sampling.iter,
#                                      alpha = alpha_set, nu = nu_set, # TODO: do we need these?
#                                      howDatSim = howDatSim,
#                                      theadaptdelta = myparams$adaptDeltaVals[8],
#                                      thetreedepth = myparams$maxTreeDepthVals[8],
#                                      thexidf = myparams$xidf_set[8],
#                                      today = todayDate,
#                                      lblabel = lblabel_set,
#                                      ublabel = ublabel_set)
# 
# outN13xi3_1iter_bayes

######################
# # multiple iterations, 1 factor combo
# 
# # choose data set (parameter combination) to work with
# oneSimSet <- simdat[simfiles[6]] # N13xi3
# runForManyRepeats(oneSimSet, 
#                   numberOfSims = n_iter,
#                   num.iters = sampling.iter,
#                   alpha = alpha_set, nu = nu_set, # TODO: do we need these?
#                   howDatSim = howDatSim,
#                   theadaptdelta = myparams$adaptDeltaVals[6],
#                   thetreedepth = myparams$maxTreeDepthVals[6],
#                   thexidf = myparams$xidf_set[6],
#                   today = todayDate,
#                   lblabel = lblabel_set,
#                   ublabel = ublabel_set)


##########################################

# multiple iterations, all factor combinations

totalTimeStart <- Sys.time()

for(j in 1:dim(myparams)[1]){
  oneSimSet <- simdat[simfiles[j]]
  runForManyRepeats(oneSimSet,
                    numberOfSims = n_iter,
                    num.iters = sampling.iter,
                    alpha = alpha_set, nu = nu_set, # TODO: do we need these?
                    howDatSim = howDatSim,
                    theadaptdelta = myparams$adaptDeltaVals[j],
                    thetreedepth = myparams$maxTreeDepthVals[j],
                    thexidf = myparams$xidf_set[j],
                    today = todayDate,
                    lblabel = lblabel_set,
                    ublabel = ublabel_set)
}

totalTime <- round(difftime(Sys.time(), totalTimeStart, units = "mins")[[1]], 3)
print(paste("Total time =", totalTime, "min"))
cat(paste("Total time =", totalTime, "min"), 
    file=paste(path, BayesCIfolder, howDatSim, "times", todayDate, ".txt", sep=""), append=T)



############################################################################
### Load output from Bayesian model
############################################################################

runDate <- "20260521" # choose date
# howDatSim <- "" # choose simulation scheme if other than Random Effects
CIfiles <- list.files(path=paste0(path, BayesCIfolder), pattern=paste0(howDatSim,"CIout",runDate)) # get file paths
CIout <- lapply(paste0(path, BayesCIfolder, CIfiles), read.csv)
# name the list elements after the N,xi parameters
Nxiparams <- c()
for (p in 1:length(CIout)) { Nxiparams <- append(Nxiparams, paste0("N",CIout[[p]][1,"N"], "xi",CIout[[p]][1,"xi"])) }
names(CIout) <- Nxiparams

############################################################################
### Check runs
############################################################################

checkBadRuns=function(CIout){
  bad1=dim(filter(CIout,parameter=="mu",min.bfmi<.3))[1] 
  bad2=dim(filter(CIout,parameter=="mu",min.bfmi<.3 | min.neff <300))[1]
  bad3=dim(filter(CIout,parameter=="mu",min.bfmi<.3 | min.neff <300 | bad.rhat > 0))[1]
  bad4=dim(filter(CIout,parameter=="mu",min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0))[1]
  bad5=dim(filter(CIout,parameter=="mu",min.bfmi<.3 | min.neff <300 | bad.rhat > 0 | num.Divergent > 0 | tree.Depth.Ex >0))[1]
  
  # print(paste("BFMI",bad1))
  # print(paste("neff",bad2-bad1))
  # print(paste("rhat",bad3-bad2))
  # print(paste("div",bad4-bad3))
  # print(paste("tree",bad5-bad4))
  # print(paste("total dropped",bad5))
  
  return(data.frame("BFMI"=bad1, "neff"=bad2-bad1, "rhat"=bad3-bad2, "div"=bad4-bad3, 
                    "tree"=bad5-bad4, "total dropped"=bad5))
}

bad.count <- bad.CIout <- final.CIout <- list()

for (p in 1:length(CIout)) {
  
  print(names(CIout)[p])
  
  # count the number of bad runs
  bad.count[[p]] <- checkBadRuns(CIout[[p]])
  
  # find the bad runs
  bad.CIout[[p]] <- filter(CIout[[p]], min.bfmi<.3 | min.neff <300 | bad.rhat > 0 
                         | num.Divergent > 0 | tree.Depth.Ex >0)
  # collect the good runs
  final.CIout[[p]] <- filter(CIout[[p]], is.na(min.bfmi) | !(min.bfmi<.3 | min.neff <300 | bad.rhat > 0
                                                          | num.Divergent > 0 | tree.Depth.Ex >0))
  # check if all runs are accounted for
  print((dim(final.CIout[[p]])[1]+dim(bad.CIout[[p]])[1]) == dim(CIout[[p]])[1])
  
}


names(bad.count) <- names(bad.CIout) <- names(final.CIout) <- Nxiparams
bad.count <- bind_rows(bad.count, .id = "parameters")
bad.count$percent.dropped <- bad.count$total.dropped/n_iter*100
bad.count

# save to csv
# write.csv(bad.count,
#           paste0(path, successmetricsfolder,
#                  "badRunsCount",howDatSim,"_",analysisMethod,"_",n_iter,"iter_",
#                  runDate,".csv"), row.names=FALSE)


############################################################################
### Success Metrics
############################################################################

# remember to use final.CIout from here on out

successMetrics <- list()

for (j in 1:length(final.CIout)) {
  p.outs <- final.CIout[[j]]
  successMetrics[[j]] <- data.frame("params" = names(final.CIout)[j],
                                    "true.mu" = p.outs$mu[1],
                                    "N" = p.outs$N[1],
                                    "true.xi" = p.outs$xi[1],
                                    "mu.bias" = mean(p.outs[p.outs$parameter=="mu", "estimate"] - p.outs$mu[1]), # avg bias
                                    "mu.covprob" = sum(p.outs[p.outs$parameter=="mu", "inInt"])/
                                                   nrow(p.outs[p.outs$parameter=="mu",]), # cp's from fitGaussGauss
                                    "xi.bias" = mean(p.outs[p.outs$parameter=="xi", "estimate"] - p.outs$xi[1]), # avg bias
                                    "xi.propbias" = mean((p.outs[p.outs$parameter=="xi", "estimate"] - 
                                                            p.outs$xi[1])/p.outs$xi[1]*100), # proportional
                                    "xi.covprob" = sum(p.outs[p.outs$parameter=="xi", "inInt"])/
                                                   nrow(p.outs[p.outs$parameter=="xi",]),
                                    "analysis.method" = analysisMethod,
                                    "sim.method" = howDatSim)
  
  # print(successMetrics[[j]])
}

successMetrics <- bind_rows(successMetrics)
successMetrics

# save to csv
# write.csv(successMetrics,
#           paste0(path, successmetricsfolder,
#                  "successMetrics",howDatSim,"_",analysisMethod,"_",n_iter,"iter_",
#                  runDate,".csv"), row.names=FALSE)
