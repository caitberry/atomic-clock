############################################################################
### 
### Dark uncertainty for atomic clock data
### Examine success metrics - plot results
### 
### Angela Folz
### April 2026
### 
############################################################################

library(readr)
library(ggplot2)

path <- "/Users/adf2/OneDrive - NIST/Documents/SpectralAnalysis/atomic-clock/"
successmetricsfolder <- "DarkUncertaintyAFST/successMetrics/"

############################################################################
### Load success metrics
############################################################################

# load from csv
successMetrics.DL <- read_csv(paste0(path, successmetricsfolder,
                                     "successMetricsRandomEffects_DL_10000iter_20260408.csv"))
successMetrics.MP <- read_csv(paste0(path, successmetricsfolder,
                                     "successMetricsRandomEffects_PM_10000iter_20260409.csv"))
successMetrics.Bayes <- read_csv(paste0(path, successmetricsfolder,
                                        "successMetricsRandomEffects_Bayes_100iter_20260428.csv"))

# choose which method to plot (ONE)
successMetrics <- successMetrics.DL; analysisMethod <- "DL"
successMetrics <- successMetrics.MP; analysisMethod <- "MP"
successMetrics <- successMetrics.Bayes; analysisMethod <- "Bayesian"

# simulation method
howDatSim <- "RandomEffects"

# number of simulation iterations
n_iter <- 100 # 10000

############################################################################
### Plot success metrics
############################################################################

figfolder <- "DarkUncertaintyAFST/figures/"

# mu bias
pdf(paste0(path, figfolder, "REsim_", analysisMethod, "_mubias_", n_iter, "iter.pdf"), width=8, height=6)

ggplot(successMetrics, aes(x=N, y=mu.bias, color=as.factor(true.xi))) +
  geom_point(size=1) +
  geom_line() +
  geom_hline(yintercept=0) +
  # ylim(c(-0.25, 0.05)) +
  theme_bw() +
  labs(color="true xi") +
  ggtitle(paste(analysisMethod, "est mu bias"))

dev.off()

# mu coverage probability
pdf(paste0(path, figfolder, "REsim_", analysisMethod, "_mucp_", n_iter, "iter.pdf"), width=8, height=6)

ggplot(successMetrics, aes(x=N, y=mu.covprob, color=as.factor(true.xi))) +
  geom_point(size=1) +
  geom_line() +
  geom_hline(yintercept=0.95) +
  # ylim(c(0.925, 0.955)) +
  theme_bw() +
  labs(color="true xi") +
  ggtitle(paste(analysisMethod, "est mu coverage probability"))

dev.off()

# xi bias
pdf(paste0(path, figfolder, "REsim_", analysisMethod, "_xibias_", n_iter, "iter.pdf"), width=8, height=6)

ggplot(successMetrics, aes(x=N, y=xi.bias, color=as.factor(true.xi))) +
  geom_point(size=1) +
  geom_line() +
  geom_hline(yintercept=0) +
  # ylim(c(-0.25, 0.05)) +
  theme_bw() +
  labs(color="true xi") +
  ggtitle(paste(analysisMethod, "est xi bias"))

dev.off()

# xi coverage probability
pdf(paste0(path, figfolder, "REsim_", analysisMethod, "_xicp_", n_iter, "iter.pdf"), width=8, height=6)

ggplot(successMetrics, aes(x=N, y=xi.covprob, color=as.factor(true.xi))) +
  geom_point(size=1) +
  geom_line() +
  geom_hline(yintercept=0.95) +
  # ylim(c(0.925, 0.955)) +
  theme_bw() +
  labs(color="true xi") +
  ggtitle(paste(analysisMethod, "est xi coverage probability"))

dev.off()


###################################
# methods on same plots

figfolder <- "DarkUncertaintyAFST/figures/compareMethods/"
successMetrics <- successMetrics.MP # don't change for this section
whichMethods <- "DLMPB"

# mu bias
pdf(paste0(path, figfolder, "REsim_", whichMethods, "_mubias_", n_iter, "iter.pdf"), width=8, height=6)

ggplot(successMetrics, aes(x=N, y=mu.bias, color=as.factor(true.xi))) +
  geom_point(size=1) +
  geom_line(aes(linetype="MP")) +
  geom_point(aes(x=successMetrics.DL$N, y=successMetrics.DL$mu.bias, color=as.factor(successMetrics.DL$true.xi))) +
  geom_line(aes(x=successMetrics.DL$N, y=successMetrics.DL$mu.bias, color=as.factor(successMetrics.DL$true.xi), 
                linetype="DL")) +
  # geom_point(aes(x=successMetrics.Bayes$N, y=successMetrics.Bayes$mu.bias, color=as.factor(successMetrics.Bayes$true.xi))) +
  # geom_line(aes(x=successMetrics.Bayes$N, y=successMetrics.Bayes$mu.bias, color=as.factor(successMetrics.Bayes$true.xi), 
  #               linetype="Bayes")) +
  geom_hline(yintercept=0) +
  # ylim(c(-0.25, 0.05)) +
  theme_bw() +
  labs(color="true xi", linetype="method") +
  ggtitle(paste("Est mu bias"))

dev.off()

# mu coverage probability
pdf(paste0(path, figfolder, "REsim_", whichMethods, "_mucp_", n_iter, "iter.pdf"), width=8, height=6)

ggplot(successMetrics, aes(x=N, y=mu.covprob, color=as.factor(true.xi))) +
  geom_point(size=1) +
  geom_line(aes(linetype="MP")) +
  geom_point(aes(x=successMetrics.DL$N, y=successMetrics.DL$mu.covprob, color=as.factor(successMetrics.DL$true.xi))) +
  geom_line(aes(x=successMetrics.DL$N, y=successMetrics.DL$mu.covprob, color=as.factor(successMetrics.DL$true.xi), 
                linetype="DL")) +
  geom_hline(yintercept=0.95) +
  ylim(c(0.925, 0.955)) +
  theme_bw() +
  labs(color="true xi", linetype="method") +
  ggtitle(paste("Est mu coverage probability"))

dev.off()

# xi bias
pdf(paste0(path, figfolder, "REsim_", whichMethods, "_xibias_", n_iter, "iter.pdf"), width=8, height=6)

ggplot(successMetrics, aes(x=N, y=xi.bias, color=as.factor(true.xi))) +
  geom_point(size=1) +
  geom_line(aes(linetype="MP")) +
  geom_point(aes(x=successMetrics.DL$N, y=successMetrics.DL$xi.bias, color=as.factor(successMetrics.DL$true.xi))) +
  geom_line(aes(x=successMetrics.DL$N, y=successMetrics.DL$xi.bias, color=as.factor(successMetrics.DL$true.xi), 
                linetype="DL")) +
  geom_hline(yintercept=0) +
  ylim(c(-0.25, 0.05)) +
  theme_bw() +
  labs(color="true xi", linetype="method") +
  ggtitle(paste("Est xi bias"))

dev.off()

# xi coverage probability
pdf(paste0(path, figfolder, "REsim_", whichMethods, "_xicp_", n_iter, "iter.pdf"), width=8, height=6)

ggplot(successMetrics, aes(x=N, y=xi.covprob, color=as.factor(true.xi))) +
  geom_point(size=1) +
  geom_line(aes(linetype="MP")) +
  geom_hline(yintercept=0.95) +
  geom_point(aes(x=successMetrics.DL$N, y=successMetrics.DL$xi.covprob, color=as.factor(successMetrics.DL$true.xi))) +
  geom_line(aes(x=successMetrics.DL$N, y=successMetrics.DL$xi.covprob, color=as.factor(successMetrics.DL$true.xi), 
                linetype="DL")) +
  ylim(c(0.925, 0.955)) +
  theme_bw() +
  labs(color="true xi", linetype="method") +
  ggtitle(paste("Est xi coverage probability"))

dev.off()

###################################
# differences between methods

successDiffs <- successMetrics.DL
successDiffs$method <- difference <- "DL-MP"
successDiffs[,c("mu.bias", "mu.covprob", "xi.bias", "xi.covprob")] <- 
  successMetrics.DL[,c("mu.bias", "mu.covprob", "xi.bias", "xi.covprob")] - 
  successMetrics.MP[,c("mu.bias", "mu.covprob", "xi.bias", "xi.covprob")]

figfolder <- "DarkUncertaintyAFST/figures/differences/"

# mu bias
pdf(paste0(path, figfolder, "REsim_", difference, "_mubias_", n_iter, "iter.pdf"), width=8, height=6)

ggplot(successDiffs, aes(x=N, y=mu.bias, color=as.factor(true.xi))) +
  geom_point(size=1) +
  geom_line() +
  geom_hline(yintercept=0) +
  # ylim(c(-0.0035, 0.0035)) +
  theme_bw() +
  labs(color="true xi") +
  ggtitle(paste(difference, "differences in est mu bias"))

dev.off()

# mu coverage probability
pdf(paste0(path, figfolder, "REsim_", difference, "_mucp_", n_iter, "iter.pdf"), width=8, height=6)

ggplot(successDiffs, aes(x=N, y=mu.covprob, color=as.factor(true.xi))) +
  geom_point(size=1) +
  geom_line() +
  geom_hline(yintercept=0) +
  ylim(c(-0.002, 0.002)) +
  theme_bw() +
  labs(color="true xi") +
  ggtitle(paste(difference, "differences in est mu coverage probability"))

dev.off()

# xi bias
pdf(paste0(path, figfolder, "REsim_", difference, "_xibias_", n_iter, "iter.pdf"), width=8, height=6)

ggplot(successDiffs, aes(x=N, y=xi.bias, color=as.factor(true.xi))) +
  geom_point(size=1) +
  geom_line() +
  geom_hline(yintercept=0) +
  ylim(c(-0.0035, 0.0035)) +
  theme_bw() +
  labs(color="true xi") +
  ggtitle(paste(difference, "differences in est xi bias"))

dev.off()

# xi coverage probability
pdf(paste0(path, figfolder, "REsim_", difference, "_xicp_", n_iter, "iter.pdf"), width=8, height=6)

ggplot(successDiffs, aes(x=N, y=xi.covprob, color=as.factor(true.xi))) +
  geom_point(size=1) +
  geom_line() +
  geom_hline(yintercept=0) +
  ylim(c(-0.002, 0.002)) +
  theme_bw() +
  labs(color="true xi") +
  ggtitle(paste(difference, "differences in est xi coverage probability"))

dev.off()
