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

path <- "/Users/adf2/OneDrive - NIST/Documents/SpectralAnalysis/atomic-clock/"
gitpath <- "Code/ComparisonAnalysis2025/DarkUncertaintyMethods/"
simdatfolder <- "DarkUncertaintyAFST/simulatedData/"
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
analyzeREff <- function(simDat, REffmethod, n_iter, showPlot){
  
  successMetrics <- list()
  
  for (j in 1:length(names(simDat))) {
    p <- names(simDat)[j]
    dat <- as.data.frame(simDat[[p]])
    true.mu <- dat$mu[1]
    true.xi <- dat$xi[1]
    Ndays <- dat$N[1]
    
    out <- list()
    ests <- list()
    
    for (i in 1:n_iter) {
      dat_iters <- dat[dat$iter==i,]
      
      ests[[i]] <- tryCatch({
        out[[i]] <- rma(yi=dat_iters$x, sei=dat_iters$u, method=REffmethod)
      
        data.frame(mu.est=out[[i]]$b, 
                   mu.cilb=out[[i]]$ci.lb, # metafor CIs all 95% by default
                   mu.ciub=out[[i]]$ci.ub, 
                   as.list(confint(out[[i]])$random["tau",]))
      }, error = function(e) {
        message(paste("Skipping index", i, "due to error:", e$message))
        return(NULL)
      })
      
    }
    ests <- bind_rows(ests)
    names(ests)[4:6] <- c("xi.est", "xi.cilb", "xi.ciub")
    row.names(ests) <- NULL
    
    if (showPlot) {
      
      # plot est means and xis
      
      p.mu <- ggplot(ests, aes(x=1:n_iter, y=mu.est)) +
        geom_point(size=1) +
        geom_errorbar(aes(ymin=mu.cilb, ymax=mu.ciub), width=0) +
        geom_hline(yintercept=true.mu) +
        geom_hline(yintercept=mean(ests$mu.est), col="orange") +
        xlab("Iteration") +
        theme_bw() +
        ggtitle(paste0(REffmethod, " est mu, N", Ndays, "xi", true.xi))
      
      p.xi <- ggplot(ests, aes(x=1:n_iter, y=xi.est)) +
        geom_point(size=1) +
        geom_errorbar(aes(ymin=xi.cilb, ymax=xi.ciub), width=0) +
        geom_hline(yintercept=true.xi) +
        geom_hline(yintercept=mean(ests$xi.est), col="orange") +
        xlab("Iteration") +
        theme_bw() +
        ggtitle(paste0(REffmethod, " est xi, N", Ndays, "xi", true.xi))
      
      print(p.mu)
      print(p.xi)
      
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
    }
    
    
    # Success metrics # TODO: save estimates, not just success metrics?????????????
    
    # coverage probabilities
    ests <- ests |> 
      mutate(mu.in = if_else(mu.cilb <= true.mu & true.mu <= mu.ciub, 1, 0))
    
    ests <- ests |> 
      mutate(xi.in = if_else(xi.cilb <= true.xi & true.xi <= xi.ciub, 1, 0))
    
    successMetrics[[j]] <- data.frame("params" = p,
                                      "true.mu" = true.mu,
                                      "N" = Ndays,
                                      "true.xi" = true.xi,
                                      "mu.bias" = mean(ests$mu.est) - true.mu, # TODO: make proportional?
                                      "mu.covprob" = sum(ests$mu.in)/n_iter,
                                      "xi.bias" = mean(ests$xi.est) - true.xi,
                                      "xi.covprob" = sum(ests$xi.in)/n_iter,
                                      "method" = REffmethod)
    
    print(successMetrics[[j]])
  }
  
  return(successMetrics)
}

############################################################################
### Analyze simulated data using random effects model
############################################################################

# choose method (ONE)
analysisMethod <- "DL"
# analysisMethod <- "PM"

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

successMetrics <- analyzeREff(simdat, analysisMethod, n_iter, showPlot=FALSE)

successMetrics <- bind_rows(successMetrics)

# save to csv
# write.csv(successMetrics,
#           paste0(path, successmetricsfolder,
#                  "successMetrics",howDatSim,"_",analysisMethod,"_",n_iter,"iter_",
#                  format(Sys.Date(), "%Y%m%d"),".csv"), row.names=FALSE)

