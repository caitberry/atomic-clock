rm(list=ls())
library(readr)
library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("MB_functions.R")
source("linearOP_ST.R")

data_path = "/Users/smt3/Documents/GitHub/atomic-clock/DarkUncertaintyAFST/simulatedData"
k_cov_factor = 1.96 

N = 13 #5, 13, 33, 100


##---LP for REM Data----------------------------------
xi = 3 #1, 3, 10
rem_data = read_csv(paste0(data_path, "/simDataRandomEffects_N", N, "xi", xi, "_10000iter_20260326.csv"))
rem_mu = rem_data$mu[1]

## Break data into a list of matricies with N rows
chunk_size = rem_data$N[1]
rem_data_trim = rem_data %>% select(x, u)
rem_data_list = split(rem_data_trim, (seq_len(nrow(rem_data_trim)) - 1) %/% chunk_size)
rem_lop_res = lapply(rem_data_list, function(input) linearOP(x = input$x, u = input$u))


coverage_rem = 0
mu_hat_rem = NULL
u_mu_hat_rem = NULL
for (i in seq_along(rem_data_list)){
  tmp = rem_lop_res[[i]]

  tmp_CI = quantile(tmp, c(0.025, 0.975))
  if((tmp_CI[1] <= rem_mu) & (tmp_CI[2] >= rem_mu)){coverage_rem = coverage_rem + 1}

  mu_hat_rem = c(mu_hat_rem, mean(tmp))
  u_mu_hat_rem = c(u_mu_hat_rem, sd(tmp))
}
coverage_rem/length(rem_data_list)

mean(mu_hat_rem)
mean(u_mu_hat_rem)


##---LP for MB Data----------------------------------
c_true = 2 #1.5, 2, 2.5

mb_data = read.table(paste0(data_path, "/simDataMulBirge_N", N, "c", c_true, "_10000iter_20260528.csv"), sep = ",", header=TRUE)
mb_mu = 0

chunk_size = N
mb_data_trim = mb_data %>% select(x, u)
mb_data_list = split(mb_data_trim, (seq_len(nrow(mb_data_trim)) - 1) %/% chunk_size)
mb_lop_res = lapply(mb_data_list, function(input) linearOP(x = input$x, u = input$u))


coverage_mb = 0
mu_hat_mb = NULL
u_mu_hat_mb = NULL
for (i in seq_along(mb_data_list)){
  tmp = mb_lop_res[[i]]

  tmp_CI = quantile(tmp, c(0.025, 0.975))
  if((tmp_CI[1] <= mb_mu) & (tmp_CI[2] >= mb_mu)){coverage_mb = coverage_mb + 1}

  mu_hat_mb = c(mu_hat_mb, mean(tmp))
  u_mu_hat_mb = c(u_mu_hat_mb, sd(tmp))
}
coverage_mb/length(mb_data_list)

mean(mu_hat_mb)
mean(u_mu_hat_mb)
