############################################################################
### 
### Dark uncertainty for atomic clock data
### Analyze simulated data assuming random effects model
### 
### Angela Folz
### April 2026
### Code adapted from Amanda Koepke
### 
############################################################################

library(readr)
library(dplyr)
library(metafor)
library(ggplot2)

path <- "/Users/adf2/OneDrive - NIST/Documents/SpectralAnalysis/atomic-clock/"
gitpath <- "Code/ComparisonAnalysis2025/DarkUncertaintyMethods/"
simdatfolder <- "DarkUncertaintyAFST/simulatedData/"
REestimatesfolder <- "DarkUncertaintyAFST/REests/"
successmetricsfolder <- "DarkUncertaintyAFST/successMetrics/"

todayDate <- format(Sys.Date(), "%Y%m%d")

############################################################################
### Load simulated data
############################################################################

# choose data to analyze (which model it was simulated under) (ONE)
# howDatSim <- "RandomEffects"
howDatSim <- "MulBirge"

simfiles <- list.files(path=paste0(path, simdatfolder), pattern=c(howDatSim,".csv")) # get file paths
simdat <- lapply(paste0(path, simdatfolder, simfiles), ifelse(howDatSim=="MulBirge", read.csv2, read.csv))
names(simdat) <- basename(simfiles) # name the list elements after the files

n_iter <- 10000 # simdat[[1]][nrow(simdat[[1]]),"iter"] # get number of simulation iterations

if(howDatSim=="MulBirge") {
  simdat[[1]]$N <- 13
  simdat[[1]]$mu <- 0
  simdat[[1]]$xi <- NA
  simdat[[1]]$iter <- rep(1:n_iter, each=simdat[[1]]$N[1])
}

############################################################################
### Functions
############################################################################

# multiple iterations
analyzeREff <- function(simDat, n_iter, REffmethod, howDatSim, today){
  
  for (j in 1:length(names(simDat))) {
    
    startTime=Sys.time()
    
    p <- names(simDat)[j]
    dat <- as.data.frame(simDat[[p]])
    true.mu <- dat$mu[1]
    true.xi <- dat$xi[1]
    Ndays <- dat$N[1]
    
    ests <- list()
    
    for (i in 1:n_iter) {
      dat_iters <- dat[dat$iter==i,]
      
     tryCatch({
        out <- rma(yi=dat_iters$x, sei=dat_iters$u, method=REffmethod)
        
        ests[[i]] <- data.frame(mu.est = out$b, 
                                mu.cilb = out$ci.lb, # metafor CIs all 95% by default
                                mu.ciub = out$ci.ub, 
                                setNames(
                                  as.list(confint(out)$random["tau",]), 
                                  c("xi.est", "xi.cilb", "xi.ciub")), # rename xi columns
                                N = Ndays, 
                                mu = true.mu, 
                                xi = true.xi)
        
        # for coverage probabilities
        ests[[i]]$mu.in = if_else(ests[[i]]$mu.cilb <= true.mu & true.mu <= ests[[i]]$mu.ciub, TRUE, FALSE)
        ests[[i]]$xi.in = if_else(ests[[i]]$xi.cilb <= true.xi & true.xi <= ests[[i]]$xi.ciub, TRUE, FALSE)
        
      }, error = function(e) {
        message(paste("Skipping index", i, "due to error:", e$message))
        return(NULL)
      })
      
      if(i %% 100==0) { # save 100 at a time
        print(paste("iteration", i))
        saveEsts <- bind_rows(ests[(i-99):i])
        write_csv(saveEsts, paste(path, REestimatesfolder, "estimates", howDatSim, REffmethod, today,
                                  "_xi",true.xi, "_N",Ndays, "_",n_iter,"iter.csv", sep=""), append=T,
                  col_names=!file.exists(paste(path, REestimatesfolder, "estimates", howDatSim, REffmethod, today,
                                               "_xi",true.xi, "_N",Ndays, "_",n_iter,"iter.csv", sep="")))
      }
    }
    endTime=Sys.time()-startTime
    print(endTime)
    cat(paste("N = ",Ndays, ", xi = ",true.xi, ", ",n_iter," iterations: run time =", endTime, "\n", sep=""), 
        file = paste(path, REestimatesfolder, howDatSim, analysisMethod, "times", today, ".txt", sep=""), append=T)
    
    return("Done!")
  }
}
   

############################################################################
### Analyze simulated data using random effects model
############################################################################

# choose method (ONE)
# analysisMethod <- "DL"
analysisMethod <- "PM"

##########################################
# # Test whether method works
# 
# # choose data set to work with
# pcombo <- simfiles[1] # simfiles[6] # choose parameter combination
# dat <- as.data.frame(simdat[[pcombo]]) # as.data.frame(simdat[pcombo])
# 
# # one iteration
# it <- 1 # choose which iteration
# dat_1iter <- dat[dat$iter==it,]
# 
# outN13xi3_1iter <- rma(yi=dat_1iter$x, sei=dat_1iter$u, method=analysisMethod)
# print(outN13xi3_1iter)
# # mu est
# outN13xi3_1iter$b # TODO: b vs beta?????????
# # mu ci
# c(outN13xi3_1iter$ci.lb, outN13xi3_1iter$ci.ub)
# # xi est
# sqrt(outN13xi3_1iter$tau2)
# # xi ci
# confint(outN13xi3_1iter)
# # confint(outN13xi3_1iter)$random["tau",]
# 
# # multiple iterations, 1 factor combo
# outN13xi3 <- analyzeREff(simdat[pcombo], analysisMethod, n_iter, showPlot=FALSE)

##########################################

# multiple iterations, all factor combinations

totalTimeStart <- Sys.time()

analyzeREff(simdat, n_iter, analysisMethod, howDatSim, todayDate)

totalTimeEnd <- Sys.time()
print(paste("Total time =", totalTimeEnd-totalTimeStart))
cat(paste("Total time =", totalTimeEnd-totalTimeStart), 
    file=paste(path, REestimatesfolder, howDatSim, analysisMethod, "times", todayDate, ".txt", sep=""), append=T)


############################################################################
### Load output from analysis
############################################################################

runDate <- "20260430" # choose date
howDatSim <- "MulBirge" # choose simulation scheme
estfiles <- list.files(path=paste0(path, REestimatesfolder), 
                       pattern=paste0("estimates",howDatSim,analysisMethod,runDate)) # get file paths
ests <- lapply(paste0(path, REestimatesfolder, estfiles), read.csv)
# name the list elements after the N,xi parameters
Nxiparams <- c()
for (p in 1:length(ests)) { Nxiparams <- append(Nxiparams, paste0("N",ests[[p]][1,"N"], "xi",ests[[p]][1,"xi"])) }
names(ests) <- Nxiparams

  
# plot est means and xis
for(j in 1:length(ests)){
  
  p.mu <- ggplot(ests[[j]][1:100,], aes(x=1:100, y=mu.est)) + # plot first 100 iterations
    geom_point(size=1) +
    geom_errorbar(aes(ymin=mu.cilb, ymax=mu.ciub), width=0) +
    geom_hline(yintercept=ests[[j]]$mu[1]) +
    geom_hline(yintercept=mean(ests[[j]]$mu.est), col="orange") +
    xlab("Iteration") +
    theme_bw() +
    ggtitle(paste0(analysisMethod, " est mu, N", ests[[j]]$N[1], "xi", ests[[j]]$xi[1]))
  print(p.mu)
  
  p.xi <- ggplot(ests[[j]][1:100,], aes(x=1:100, y=xi.est)) + # plot first 100 iterations
    geom_point(size=1) +
    geom_errorbar(aes(ymin=xi.cilb, ymax=xi.ciub), width=0) +
    geom_hline(yintercept=ests[[j]]$xi[1]) +
    geom_hline(yintercept=mean(ests[[j]]$xi.est), col="orange") +
    xlab("Iteration") +
    theme_bw() +
    ggtitle(paste0(analysisMethod, " est xi, N", ests[[j]]$N[1], "xi", ests[[j]]$xi[1]))
  print(p.xi)
}
  
  # histograms
  # h.mu <- ggplot(ests, aes(x=mu.est)) +
  #   geom_histogram() +
  #   geom_vline(xintercept=true.mu) +
  #   theme_bw() +
  #   ggtitle(paste0(REffmethod, " est mu, N", Ndays, "xi", true.xi))
  # 
  # h.xi <- ggplot(ests, aes(x=xi.est)) +
  #   geom_histogram() +
  #   geom_vline(xintercept=true.xi) +
  #   theme_bw() +
  #   ggtitle(paste0(REffmethod, " est xi, N", Ndays, "xi", true.xi))
  # 
  # print(h.mu)
  # print(h.xi)


############################################################################
### Success Metrics
############################################################################

successMetrics <- list()

for (j in 1:length(ests)) {
  p.outs <- ests[[j]]
  successMetrics[[j]] <- data.frame("params" = names(ests)[j],
                            "true.mu" = p.outs$mu[1],
                            "N" = p.outs$N[1],
                            "true.xi" = p.outs$xi[1],
                            "mu.bias" = mean(p.outs$mu.est - p.outs$mu[1]), # avg bias
                            # "mu.propbias" = mean((p.outs$mu.est - p.outs$mu[1])/p.outs$mu[1]*100), # can't if mu=0
                            "mu.covprob" = sum(p.outs$mu.in)/n_iter, # cp's from fitGausGaus
                            "xi.bias" = mean(p.outs$xi.est - p.outs$xi[1]), # avg bias
                            "xi.propbias" = mean((p.outs$xi.est - p.outs$xi[1])/p.outs$xi[1]*100), # proportional
                            "xi.covprob" = sum(p.outs$xi.in)/n_iter,
                            "analysis.method" = analysisMethod,
                            "sim.method" = howDatSim)
}
successMetrics <- bind_rows(successMetrics)
successMetrics

# save to csv
# write.csv(successMetrics,
#           paste0(path, successmetricsfolder,
#                  "successMetrics",howDatSim,"_",analysisMethod,"_",n_iter,"iter_",
#                  runDate,".csv"), row.names=FALSE)

