

# Al+/Sr -0.2 ± 1.9 -0.3 ± 2.0 1.0 (9 days) (10^-18 scale, i'm going to work in -16)
9 days
mu around: -200
tau around: 1
daily uncertainties: ??


# Al+/Yb 0.19 ± 3.0 0.2 ± 3.0 0.9 (9 days)
9 days
mu around: .2
tau around: 1
daily uncertainties: ??
  
# Yb/Sr -0.5 ± 3.0 -0.7 ± 3.0 2.9 (13 days)
13 days
mu around: -.5
tau around: 3
daily uncertainties:
  
# --- 1. Load Data ---
df_AlSr <- read_csv("Data/ClockComparison2025/BayesianAnalysisData/ErYb_AlSr_data.csv")
df_AlYb <- read_csv("Data/ClockComparison2025/BayesianAnalysisData/ErYb_AlYb_data.csv")
df_YbSr <- read_csv("Data/ClockComparison2025/BayesianAnalysisData/ErYb_YbSr_data.csv")


data=df_AlSr
mu1=-200
tau1=1
n1=9
sigma1=runif(n1,min=min(data$statistical_unc),max = max(data$statistical_unc))

sim1=rnorm(n1,mu1,tau1)


