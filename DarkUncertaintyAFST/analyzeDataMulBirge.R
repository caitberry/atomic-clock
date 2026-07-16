############################################################################
### 
### Dark uncertainty for atomic clock data
### Analyse data with MB regardless of how it was simulated 
### 
### Suzanne Thornton
### July 2026
### Code adapted from Amanda Koepke and Angela Folz
### Code optimization assistance from Gemini
### 
############################################################################

rm(list=ls())
library(data.table)
library(tidyverse)
library(gridExtra)

path = "/Users/smt3/Documents/GitHub/atomic-clock/"
simdatfolder = "DarkUncertaintyAFST/simulatedData/" #added to gitignore
figfolder = "DarkUncertaintyAFST/figures/" #added to gitignore

##Supporter functions
source(paste0(path, "DarkUncertaintyAFST/", "MB_functions.R"))

##BACON data
ratiolab = "YbSr"
ratiodf = read_csv(paste0(path, "Data/ClockComparison2025/BayesianAnalysisData/ErYb_",ratiolab,"_data.csv"))
bacon_measurements = ratiodf$offset
bacon_uncertainties = ratiodf$statistical_unc

N_cov = 1E4
k_cov_factor = 1.96 
N = length(bacon_uncertainties)

##---------------------------------------------------------------
##---Analyse MB Data for N=13 for table
##---------------------------------------------------------------N = length(bacon_measurements)
birg_constant = 2
mu_MB = 0


## Simulated data (Vectorized)
uncertainties = runif(N, min(bacon_uncertainties), max(bacon_uncertainties)) 
measurements = rnorm(N, mean = mu_MB, sd = uncertainties * birg_constant)
sim_dat = data.frame(day = 1:N, x = measurements, u = uncertainties)
summary(sim_dat)

mu_hat_b = weighted_mean(sim_dat$x, sim_dat$u)
u_mu_hat_b = u_MB(sim_dat$x, sim_dat$u)
u_mu_hat_b_corrected = u_MB(sim_dat$x, sim_dat$u, correction = TRUE)


# plot real BACON2 data
ggplot(ratiodf, aes(x=as.character(date), y=offset)) +
  geom_point(size=1) +
  geom_errorbar(aes(ymin=offset-statistical_unc, ymax=offset+statistical_unc), width=0) +
  ylim(c(-110,-90)) +
  theme_bw() +
  ggtitle("BACON2 data")

# plot simulated data set
# plot includes mu_hat with uncertainty bars showing k*u, where k=1 and true mu (red)
p1 <- ggplot(sim_dat, aes(x=day, y=x)) +
  geom_point(size=1) +
  geom_errorbar(aes(ymin=x-u, ymax=x+u), width=0) +
  geom_hline(aes(yintercept = mu_MB, color = "true mu", linetype = "true mu")) + 
  geom_hline(aes(yintercept = mu_hat_b, color = "weighted mean", linetype = "weighted mean")) + 
  geom_hline(aes(yintercept = mu_hat_b + u_mu_hat_b, color = "unadjusted Birge bounds", linetype = "unadjusted Birge bounds")) + 
  geom_hline(aes(yintercept = mu_hat_b - u_mu_hat_b, color = "unadjusted Birge bounds", linetype = "unadjusted Birge bounds")) +
  geom_hline(aes(yintercept = mu_hat_b + u_mu_hat_b_corrected, color = "adjusted Birge bounds", linetype = "adjusted Birge bounds")) + 
  geom_hline(aes(yintercept = mu_hat_b - u_mu_hat_b_corrected, color = "adjusted Birge bounds", linetype = "adjusted Birge bounds")) + 
  scale_color_manual(name = "Legend",
              breaks = c("true mu", 
               "weighted mean", 
               "unadjusted Birge bounds", 
               "adjusted Birge bounds"),
              values = c("true mu" = "red", 
                "weighted mean" = "purple", 
                "unadjusted Birge bounds" = "purple", 
                "adjusted Birge bounds" = "purple"),
              labels = c("true mu" = expression(mu), 
               "weighted mean" = expression(hat(mu)[WM]),
               "unadjusted Birge bounds" = "unadjusted Birge bounds",
               "adjusted Birge bounds" = "adjusted Birge bounds")) +
  scale_linetype_manual(name = "Legend",
              breaks = c("true mu", 
               "weighted mean", 
               "unadjusted Birge bounds", 
               "adjusted Birge bounds"),
              values = c("true mu" = 1, 
                "weighted mean" = 1, 
                "unadjusted Birge bounds" = 2, 
                "adjusted Birge bounds" = 3),
              labels = c("true mu" = expression(mu), 
               "weighted mean" = expression(hat(mu)[WM]),
               "unadjusted Birge bounds" = "unadjusted Birge bounds",
               "adjusted Birge bounds" = "adjusted Birge bounds")) +
  theme_bw() +
  ylab(expression(x[i] %+-% u(x[i]))) + xlab("Day") + 
  ggtitle(paste("Simulated data, N =", N, ", c =", birg_constant))

#ggsave("DarkUncertaintyAFST/mul_birge_plot.png", plot = p1, width = 6, height = 4, units = "in", dpi = 300)


##coverage and bias
uncertainties_cov = runif(N*N_cov, min(bacon_uncertainties), max(bacon_uncertainties))
measurements_cov = rnorm(N * N_cov, mean = mu_MB, sd = uncertainties_cov * 2)

dt_single = data.table(
  run_id = rep(1:N_cov, each = N),
  day = rep(1:N, N_cov),
  u = uncertainties_cov,
  x = measurements_cov
)

birge_res_single = dt_single[, birge_wrapper(as.matrix(.SD)), by = run_id, .SDcols = c("x", "u")]

single_stats = birge_res_single[, .(
  cov_unadj = mean((mean_birge - k_cov_factor * u_birge <= mu_MB) & 
                     (mean_birge + k_cov_factor * u_birge >= mu_MB)),
  cov_adj = mean((mean_birge - k_cov_factor * u_birge_corrected <= mu_MB) & 
                   (mean_birge + k_cov_factor * u_birge_corrected >= mu_MB)),
  mean_mu_hat = mean(mean_birge),
  mean_u_unadj = mean(u_birge),
  mean_u_adj = mean(u_birge_corrected)
)]

print(single_stats)

##---------------------------------------------------------------
##---Analyse REM Data for N=13 for table
##---------------------------------------------------------------
## TODO: make more efficient using data.table as below 
rem_data = read_csv(paste0(path, "DarkUncertaintyAFST/simulatedData/simDataRandomEffects_N13xi3_10000iter_20260326.csv"))
mu_REM = rem_data$mu[1]

## Break data into a list of matricies with 13 rows
chunk_size = rem_data$N[1]
rem_data_list = split(rem_data, (seq_len(nrow(rem_data)) - 1) %/% chunk_size)
rem_mtx_list = lapply(rem_data_list, as.matrix)
rem_birge_res = lapply(rem_mtx_list, birge_wrapper) ##Note input must be a matrix

coverage = 0
coverage_corrected = 0
mu_hat_un = NULL
mu_hat_adj = NULL
u_mu_hat_un = NULL 
u_mu_hat_adj = NULL
for (i in seq_along(rem_data_list)){
  tmp = rem_birge_res[[i]]
  if((tmp$mean_birge - k_cov_factor*tmp$u_birge <= mu_REM) & (tmp$mean_birge + k_cov_factor*tmp$u_birge >= mu_REM)){ coverage = coverage + 1}
  if((tmp$mean_birge - k_cov_factor*tmp$u_birge_corrected <= mu_REM) & (tmp$mean_birge + k_cov_factor*tmp$u_birge_corrected >= mu_REM)){ coverage_corrected = coverage_corrected + 1}

  mu_hat_un = c(mu_hat_un, tmp$mean_birge)
  mu_hat_adj = c(mu_hat_adj, tmp$mean_birge)

  u_mu_hat_un = c(u_mu_hat_un, tmp$u_birge)
  u_mu_hat_adj = c(u_mu_hat_adj, tmp$u_birge_corrected)
}
coverage/length(rem_data_list)
coverage_corrected/length(rem_data_list)

mean(mu_hat_un)
mean(u_mu_hat_un)
mean(u_mu_hat_adj)

# ---------------------------------------------------------
# --- Read CSV files and produce Bias/Coverage Plots
# --- for MB analysis of MB data and REM data
# ---------------------------------------------------------
N_new = c(5, 13, 33, 100)
#mu = 0
#N_cov = 1E4 


##---MB Data--------------------------------------------------
true_c = c(1.5, 2, 3) ##Note: these must match with simulated data 
mb_data_names = "simDataMulBirge_"
mb_mu = 0

mb_data_mb_analysis = replicate(length(N_new)*length(true_c),
                        list(N = NULL, c = NULL, cov_un = NULL, cov = NULL, bias_both = NULL),
                        simplify = FALSE)
counter = 1

for (n in N_new){
  for (c in true_c){
    mb_data = read_csv(paste0(path, simdatfolder, mb_data_names, "N", n, "c", c, "_10000iter_20260528.csv"))
    
    # Extract only one row per simulation run (since day goes 1 to n)
    mb_data_summary <- mb_data %>% filter(day == 1)

    # Vectorized coverage and bias calculation (No loop required!)
    coverage = sum(mb_data_summary$lb <= mb_mu & mb_data_summary$ub >= mb_mu)
    coverage_corrected = sum(mb_data_summary$lb_corrected <= mb_mu & mb_data_summary$ub_corrected >= mb_mu)
    mu_hat = (mb_data_summary$lb + mb_data_summary$ub) / 2
    
    # Store results
    mb_data_mb_analysis[[counter]] = list(
      N = n,
      c = c,
      cov_un = coverage / nrow(mb_data_summary),
      cov_adj = coverage_corrected / nrow(mb_data_summary),
      bias = mean(mu_hat)
    )
    counter = counter + 1
  }
}

##---REM Data----------------------------------
true_xi = c(1, 2, 3) ##Note: these must match with simulated data
rem_data_names = "simDataRandomEffects_"

rem_data_mb_analysis = replicate(length(N_new)*length(true_xi),
                        list(N = NULL, xi = NULL, cov_un = NULL, cov = NULL, bias_both = NULL),
                        simplify = FALSE)
counter = 1

for (n in N_new){
  for (xi in true_xi){
    folder_dir   = paste0(path, simdatfolder)
    file_pattern = paste0("^", rem_data_names, "N", n, "xi", xi, ".*\\.csv$")
    file_path = list.files(path = folder_dir, pattern = file_pattern, full.names = TRUE)
    
    # Load directly as data.table
    rem_dt = fread(file_path[1])
    mu_REM = rem_dt$mu[1]
    chunk_size = rem_dt$N[1]

    # Old ineffecient code:
    # rem_data_list = split(rem_data, (seq_len(nrow(rem_data)) - 1) %/% chunk_size)
    # rem_mtx_list = lapply(rem_data_list, as.matrix)
    # rem_birge_res = lapply(rem_mtx_list, birge_wrapper) ##Note input must be a matrix

    # New efficient code:
    # 1. Instantly assign run IDs
    # creates one new column (run_id) which is filled by calculating a group number for every individual row in the data table 
    # .N is total number of rows in the table 
    # - 1 to shift to 0-indexed sequence for %/% which is integer division (rounds down)
    # + 1 to shift back to 1-index
    rem_dt[, run_id := (seq_len(.N) - 1) %/% chunk_size + 1]
    
    # 2. Vectorized calculation of Birge parameters by group
    rem_birge_res = rem_dt[, birge_wrapper(as.matrix(.SD)), by = run_id, .SDcols = c("x", "u")]

    # 3. Vectorized coverage and bias 
    cov_un = rem_birge_res[, mean((mean_birge - k_cov_factor * u_birge <= mu_REM) & 
                                  (mean_birge + k_cov_factor * u_birge >= mu_REM))]
    cov_adj = rem_birge_res[, mean((mean_birge - k_cov_factor * u_birge_corrected <= mu_REM) & 
                                   (mean_birge + k_cov_factor * u_birge_corrected >= mu_REM))]
    mean_bias = mean(rem_birge_res$mean_birge)

    rem_data_mb_analysis[[counter]] = list(
      N = n,
      xi = xi,
      cov_un = cov_un,
      cov_adj = cov_adj,
      bias = mean_bias
    )

    counter = counter + 1
 
  }
}


##--- Plots -------------------------------------------
# MB data on left, REM data on right
# mu coverage top, mu bias bottom
# colors for true c and true xi 
# dashed/solid lines for adjusted/unadjusted MB (ie analysis method)
# coverage prob (y-axis)
# bias (y axis)
# N = c(5, 13, 33, 100) (x-axis)

## Create data frames of info

## Old code
# mb_data_df = data.frame(N = 0, c = 0, cov_un = 0, cov_adj = 0, bias = 0)
# rem_data_df = data.frame(N = 0, xi = 0, cov_un = 0, cov_adj = 0, bias = 0)

# for (i in 1:12) {
#     mb_data_df[i, ] = c(mb_data_mb_analysis[[i]]$N,
#                         mb_data_mb_analysis[[i]]$c,
#                         mb_data_mb_analysis[[i]]$cov_un,
#                         mb_data_mb_analysis[[i]]$cov,
#                         mb_data_mb_analysis[[i]]$bias_both
#                         )
#     rem_data_df[i, ] = c(rem_data_mb_analysis[[i]]$N,
#                         rem_data_mb_analysis[[i]]$xi,
#                         rem_data_mb_analysis[[i]]$cov_un,
#                         rem_data_mb_analysis[[i]]$cov,
#                         rem_data_mb_analysis[[i]]$bias_both
#                         )
# }

## New code using data.table
mb_data_df  <- data.table::rbindlist(mb_data_mb_analysis)
rem_data_df <- data.table::rbindlist(rem_data_mb_analysis)

## Convert data frames to long format
mb_all_long <- mb_data_df %>%
  pivot_longer(
    cols = c(cov_un, cov_adj),
    names_to = "Method",
    values_to = "cov_prob"
  )
rem_all_long <- rem_data_df %>%
  pivot_longer(
    cols = c(cov_un, cov_adj),
    names_to = "Method",
    values_to = "cov_prob"
  )


## Generate and save plots
mb_cov <- ggplot(mb_all_long, aes(x=N, y=cov_prob, color=as.factor(c), linetype=as.factor(Method))) +
  geom_point(size=1) +
  geom_line() +
  geom_hline(yintercept=0.95) +
  theme_bw() +
  labs(color="True c", linetype="Method") +
  ggtitle("MB data, MB est mu coverage probability")

mb_bias <- ggplot(mb_data_df, aes(x=N, y=bias, color=as.factor(c))) +
  geom_point(size=1) +
  geom_line() +
  theme_bw() +
  labs(color="True c") +
  ggtitle("MB data, MB est mu bias")

rem_cov <- ggplot(rem_all_long, aes(x=N, y=cov_prob, color=as.factor(xi), linetype=as.factor(Method))) +
  geom_point(size=1) +
  geom_line() +
  geom_hline(yintercept=0.95) +
  theme_bw() +
  labs(color="True xi", linetype="Method") +
  ggtitle("REM data, MB est mu coverage probability")

rem_bias <- ggplot(rem_data_df, aes(x=N, y=bias, color=as.factor(xi))) +
  geom_point(size=1) +
  geom_line() +
  theme_bw() +
  labs(color="True xi") +
  ggtitle("REM data, MB est mu bias")



#pdf(paste0(path, figfolder, "sim_MB_mu_cov_bias_1E4iter.pdf"), width=8, height=6)
grid.arrange(mb_cov, rem_cov, mb_bias, rem_bias, nrow=2)
#dev.off()