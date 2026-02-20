############################################################################
### 
### Dark uncertainty for atomic clock data
### Simulate data assuming random effects model
### 
### Angela Folz
### January 2026
### Code adapted from Amanda Koepke
### 
############################################################################

library(readr)
library(dplyr)
library(ggplot2)
library(metafor)

path <- "/Users/adf2/OneDrive - NIST/Documents/SpectralAnalysis/atomic-clock/"
simdatfolder <- "DarkUncertaintyAFST/simulatedData/"

############################################################################
### Set parameters based on data
############################################################################

# ratiolab="AlSr"
# ratiolab="AlYb"
ratiolab="YbSr"

# read in means and uncertainties for multiple days of BACON2 ratio measurements
ratiodf <- read_csv(paste0(path, "Data/ClockComparison2025/BayesianAnalysisData/ErYb_",ratiolab,"_data.csv"))

measurements <- ratiodf$offset
uncertainties <- ratiodf$statistical_unc

N <- length(measurements)
mu <- 0 # centered
      # mean(measurements) # observed mean over multiple days from BACON2 ratio measurements

N_set <- c(N,N+20,100) # number of days; vary N to reflect realistic and ideal data sizes # TODO: j=1,...,9,25?
xi_set <- c(3,10) # between-day variability # TODO: use values we expect to see from the xi reported in BACON2
                                            # and one way bigger than what we see, one close to 0, one scale of YbSr
allparams <- expand.grid(N_set=N_set, xi_set=xi_set) # want every combo of N and xi

lb <- min(uncertainties)
ub <- max(uncertainties)

n_iter <- 100 #00 ### make bigger
baseSeed <- 1000 ### make bigger if increase n_iter

howDatSim <- "RandomEffects"

# test
N <- N_set[1]
xi <- xi_set[1]


############################################################################
### Simulate data from random effects model
############################################################################

simulateDataRandEff <- function(N, mu, xi, lb, ub, seed){
  
  set.seed(seed)
  lambda <- rnorm(N, 0, xi) # day effects
  
  sigma <- runif(N, lb, ub) # simulate "known" uncertainties for each day from a uniform, 
                            # with the lower and upper bound set to be what was observed 
                            # in the BACON2 data for a particular ratio
  epsilon <- rnorm(N, 0, sigma)
  
  x <- mu + lambda + epsilon # random effects model
  
  dat <- data.frame(Day=as.factor(1:N), x=x, u=sigma, 
                    N=N, mu=mu, xi=xi, lb=lb, ub=ub, seed=seed)
  return(dat)
}


############################################################################
### Run
############################################################################

# run once for one N, xi combo
oneSimDat <- simulateDataRandEff(N=N, mu=mu, xi=xi, lb=lb, ub=ub, seed=baseSeed)


# run for all N, xi combinations, n_iter times
simDat <- list()
for (p in 1:dim(allparams)[1]) {
  simDatSet <- list()
  
    for (i in 1:n_iter) {
      simDatSet[[i]] <- simulateDataRandEff(N=allparams$N_set[p], mu=mu, xi=allparams$xi_set[p],
                                             lb=lb, ub=ub, seed=baseSeed*p+i) 
                                    # TODO: double check seeds don't repeat with higher n_iter!!!!!!!!!!!
      simDatSet[[i]]$iter <- i
    }
  
  simDat[[p]] <- bind_rows(simDatSet)
  names(simDat)[[p]] <- paste0("N", allparams$N_set[p], "xi", allparams$xi_set[p])
}

# save to csv
# for (i in 1:length(simDat)) {
#   write.csv(simDat[[i]], paste0(path, simdatfolder,
#                                 "simData",howDatSim,"_",names(simDat)[i],"_",n_iter,"iter_",
#                                 format(Sys.Date(), "%Y%m%d"),".csv"))
# }


############################################################################
### Visualize
############################################################################

# Plot consensus mean with uncertainty bars showing k*u, where k=1

# plot real BACON2 data
ggplot(ratiodf, aes(x=date, y=offset)) +
  geom_point(size=1) +
  geom_errorbar(aes(ymin=offset-statistical_unc, ymax=offset+statistical_unc), width=0) +
  ylim(c(-110,-90)) +
  theme_bw() +
  ggtitle("BACON2 data")


# plot for one N, xi combo (that was run separately above)
ggplot(oneSimDat, aes(x=Day, y=x)) +
  geom_point(size=1) +
  geom_errorbar(aes(ymin=x-u, ymax=x+u), width=0) +
  ylim(c(-6,6)) +
  # ylim(c(-110,-90)) + # if use mu=mean(BACON2 data)
  theme_bw() +
  ggtitle(paste("Simulated data, N =", N, ", xi =", xi))


# plot for one iteration of one N, xi combo (from big sim data set)

names(simDat) # see which parameter combinations are available

pcombo <- "N13xi3" # choose parameter combination
dat <- as.data.frame(simDat[pcombo])
names(dat) <- c("Day", "x", "u", "N", "mu", "xi", "lb","ub","seed", "iter") # fix names
it <- 1 # choose which iteration

ggplot(dat[dat$iter==it,], aes(x=Day, y=x)) +
  geom_point(size=1) +
  geom_errorbar(aes(ymin=x-u, ymax=x+u), width=0) +
  ylim(c(-10,10)) +
  # ylim(c(-130,-70)) + # if use mu=mean(BACON2 data)
  theme_bw() +
  ggtitle(paste0("Simulated data, ", pcombo, " #", it))


# plot for multiple iterations of  all parameter combinations

iters <- c(1, 19, 72) # choose which iterations

for (p in names(simDat)) {
  dat <- as.data.frame(simDat[p])
  names(dat) <- c("Day", "x", "u", "N", "mu", "xi", "lb","ub","seed", "iter") # fix names
  
  for (it in iters) {
    g <- ggplot(dat[dat$iter==it,], aes(x=Day, y=x)) +
      geom_point(size=1) +
      geom_errorbar(aes(ymin=x-u, ymax=x+u), width=0) +
      ylim(c(-30,30)) +
      # ylim(c(-130,-70)) + # if use mu=mean(BACON2 data)
      theme_bw() +
      ggtitle(paste0("Simulated data, ", p, " #", it))
    print(g)
  }
}


############################################################################
### Load simulated data
############################################################################

# load from csv
# datN13xi3 <- read_csv(paste0(path, simdatfolder,
#                              "simDataRandomEffects_N13xi3_100iter_20260130.csv"))
# datN13xi10 <- read_csv(paste0(path, simdatfolder,
#                               "simDataRandomEffects_N13xi10_100iter_20260130.csv"))
# datN33xi3 <- read_csv(paste0(path, simdatfolder,
#                              "simDataRandomEffects_N33xi3_100iter_20260130.csv"))
# datN33xi10 <- read_csv(paste0(path, simdatfolder,
#                               "simDataRandomEffects_N33xi10_100iter_20260130.csv"))
# datN100xi3 <- read_csv(paste0(path, simdatfolder,
#                               "simDataRandomEffects_N100xi3_100iter_20260130.csv"))
# datN100xi10 <- read_csv(paste0(path, simdatfolder,
#                                "simDataRandomEffects_N100xi10_100iter_20260130.csv"))
# 
# # choose data set to work with
# dat <- datN13xi3

# OR load from stuff generated above
pcombo <- "N13xi3" # choose parameter combination
dat <- as.data.frame(simDat[pcombo])
names(dat) <- c("Day", "x", "u", "N", "mu", "xi", "lb","ub","seed", "iter") # fix names


############################################################################
### Analyze simulated data using random effects model
############################################################################


# one iteration
it <- 1 # choose which iteration
dat_1iter <- dat[dat$iter==it,]

outN13xi3_1iter <- rma(yi=dat_1iter$x, sei=dat_1iter$u, method="DL")
print(outN13xi3_1iter)
mu.est_outN13xi3_1iter <- outN13xi3_1iter$b # TODO: b vs beta?????????
mu.ci_outN13xi3_1iter <- c(outN13xi3_1iter$ci.lb, outN13xi3_1iter$ci.ub)
xi.est_outN13xi3_1iter <- sqrt(outN13xi3_1iter$tau2)
confint(outN13xi3_1iter)
confint(outN13xi3_1iter)$random["tau",]


# multiple iterations

showPlot <- TRUE # plot the estimates and CIs?

successMetrics <- list()

for (j in 1:length(names(simDat))) {
  p <- names(simDat)[j]
  dat <- as.data.frame(simDat[p])
  names(dat) <- c("Day", "x", "u", "N", "mu", "xi", "lb","ub","seed", "iter") # fix names
  true.mu <- dat$mu[1]
  true.xi <- dat$xi[1]
  Ndays <- dat$N[1]
  
  # iters <- c(1, 19, 72) # choose subset of iterations
  DLout <- list()
  DLests <- list()
  
  for (i in 1:n_iter) { # length(iters)) {
    dat_iters <- dat[dat$iter==i,] # iters[i],]
    DLout[[i]] <- rma(yi=dat_iters$x, sei=dat_iters$u, method="DL")
    DLests[[i]] <- data.frame(mu.est=DLout[[i]]$b, 
                            mu.cilb=DLout[[i]]$ci.lb, # metafor CIs all 95% by default
                            mu.ciub=DLout[[i]]$ci.ub, 
                            as.list(confint(DLout[[i]])$random["tau",]))
  }
  DLests <- bind_rows(DLests)
  names(DLests)[4:6] <- c("xi.est", "xi.cilb", "xi.ciub")
  row.names(DLests) <- NULL
  
  if (showPlot) {
  
    # plot est means and xis
    
    p.mu <- ggplot(DLests, aes(x=1:n_iter, y=mu.est)) +
      geom_point(size=1) +
      geom_errorbar(aes(ymin=mu.cilb, ymax=mu.ciub), width=0) +
      geom_hline(yintercept=true.mu) +
      geom_hline(yintercept=mean(DLests$mu.est), col="orange") +
      xlab("Iteration") +
      theme_bw() +
      ggtitle(paste0("DL est mu, N", Ndays, "xi", true.xi))
    
    p.xi <- ggplot(DLests, aes(x=1:n_iter, y=xi.est)) +
      geom_point(size=1) +
      geom_errorbar(aes(ymin=xi.cilb, ymax=xi.ciub), width=0) +
      geom_hline(yintercept=true.xi) +
      geom_hline(yintercept=mean(DLests$xi.est), col="orange") +
      xlab("Iteration") +
      theme_bw() +
      ggtitle(paste0("DL est xi, N", Ndays, "xi", true.xi))
    
    print(p.mu)
    print(p.xi)
    
    # histograms
    # h.mu <- ggplot(DLests, aes(x=mu.est)) +
    #   geom_histogram() +
    #   geom_vline(xintercept=true.mu) +
    #   theme_bw() +
    #   ggtitle(paste0("DL est mu, N", Ndays, "xi", true.xi))
    # 
    # h.xi <- ggplot(DLests, aes(x=xi.est)) +
    #   geom_histogram() +
    #   geom_vline(xintercept=true.xi) +
    #   theme_bw() +
    #   ggtitle(paste0("DL est xi, N", Ndays, "xi", true.xi))
    # 
    # print(h.mu)
    # print(h.xi)
  }


  # Success metrics
  
  # coverage probabilities
  DLests <- DLests |> 
    mutate(mu.in = if_else(mu.cilb <= true.mu & true.mu <= mu.ciub, 1, 0))
  
  DLests <- DLests |> 
    mutate(xi.in = if_else(xi.cilb <= true.xi & true.xi <= xi.ciub, 1, 0))
  
  successMetrics[[j]] <- data.frame("params" = p, 
                                "true.mu" = true.mu,
                                "N" = Ndays,
                                "true.xi" = true.xi,
                                "mu.bias" = mean(DLests$mu.est) - true.mu,
                                "mu.covprob" = sum(DLests$mu.in)/n_iter,
                                "xi.bias" = mean(DLests$xi.est) - true.xi,
                                "xi.covprob" = sum(DLests$xi.in)/n_iter)
  
  print(p) # parameter combination
  
  cat("mu bias =", mean(DLests$mu.est) - true.mu, "\n")
  cat("mu coverage probability =", sum(DLests$mu.in)/n_iter, "\n")
  
  cat("xi bias =", mean(DLests$xi.est) - true.xi, "\n")
  cat("xi coverage probability =", sum(DLests$xi.in)/n_iter, "\n")
}

successMetrics <- bind_rows(successMetrics)


# plot success metrics

figfolder <- "DarkUncertaintyAFST/figures/"

# mu bias
pdf(paste0(path, figfolder, "REsim_DL_mubias_100iter.pdf"), width=8, height=6)

ggplot(successMetrics, aes(x=N, y=mu.bias, color=as.factor(true.xi))) +
  geom_point(size=1) +
  geom_line() +
  theme_bw() +
  labs(color="true xi") +
  ggtitle("DL est mu bias")

dev.off()

# mu coverage probability
pdf(paste0(path, figfolder, "REsim_DL_mucp_100iter.pdf"), width=8, height=6)

ggplot(successMetrics, aes(x=N, y=mu.covprob, color=as.factor(true.xi))) +
  geom_point(size=1) +
  geom_line() +
  geom_hline(yintercept=0.95) +
  theme_bw() +
  labs(color="true xi") +
  ggtitle("DL est mu coverage probability")

dev.off()

# xi bias
pdf(paste0(path, figfolder, "REsim_DL_xibias_100iter.pdf"), width=8, height=6)

ggplot(successMetrics, aes(x=N, y=xi.bias, color=as.factor(true.xi))) +
  geom_point(size=1) +
  geom_line() +
  theme_bw() +
  labs(color="true xi") +
  ggtitle("DL est xi bias")

dev.off()

# xi coverage probability
pdf(paste0(path, figfolder, "REsim_DL_xicp_100iter.pdf"), width=8, height=6)

ggplot(successMetrics, aes(x=N, y=xi.covprob, color=as.factor(true.xi))) +
  geom_point(size=1) +
  geom_line() +
  geom_hline(yintercept=0.95) +
  theme_bw() +
  labs(color="true xi") +
  ggtitle("DL est xi coverage probability")

dev.off()


