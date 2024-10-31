# Title: Simulation Studies
# File name: Simulation.R
# Description: This file contains the code required to run simulations of
#              known noise processes in order to compare spectral-based 
#              AVAR and OAVAR


# libraries #############################################

mtse <- modules::use("Functions.R")


# Math ##################################################
# for a given frequency f, tau, and spectral estimate S.hat we can get an allan variance
# estimate sigma.hat(tau) using:

# sigma.hat(tau) = int_{0}^{1/2} G_tau(f)*S.hat(f)df ~ delta.f * G.vec %*% S.hat (or some equivalent quadrature)

# where 
# G.vec is G_tau evaluated at all fourier frequencies f_0, f_1, ..., f_N/2
# S.hat is the MTSE at the vector of fourier frequencies
# delta.f is the 

# Example ###############################################

#### (->) inputs to simulation

# placement of gaps

# set parameters of MTSE: W, K

# tau values of interest


#### needed values and vectors/matrices for saving results

# length of time series (without gaps)

# length of time series (with gaps)

# 




