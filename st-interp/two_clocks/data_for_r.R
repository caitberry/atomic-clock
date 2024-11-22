remove(list = ls())
library(tidyverse)
library(openxlsx)
library(imputeTS)
library(pracma)
options(digits = 22)

ErYb_data <- read.table("ErYb_data.csv", header = TRUE, sep = ",",
                        numerals = "no.loss")
#head(ErYb_data)
ErYb_data <- ErYb_data[-1, ]
#head(ErYb_data)
comb <- read.table("Common_mjd_comb.csv", header = TRUE,
                   sep=",", colClasses = c("numeric"), numerals = "no.loss")
sr_shift <- read.table("Sr_shift_data.csv", header = TRUE,
                       sep = ",", colClasses = c("numeric", "numeric"),
                       numerals = "no.loss")
yb_shift <- read.table("Yb_shift_data.csv", header = TRUE,
                       sep = ",",
                       colClasses = c("numeric", "numeric", "character"),
                       numerals = "no.loss")

head(comb)
head(sr_shift)
print(sr_shift[1, 2], digits = 22); #as.double(sr_shift[1,2], length=30)
head(yb_shift)

################################################################################
# Select only common values of MJD for analysis
################################################################################

## LB = comb$MJD[1] = sr_shift$MJD[325] = yb_shift[140]
## UB = comb$MJD[22064] = sr_shift$MJD[3126] = yb_shift$MJD[25396]
common_mjd <- comb$MJD[1:22063]
sr_shift <- sr_shift[1:3126, ]
#yb_shift <- yb_shift[140:25397, ]

common_mjd[1]; sr_shift$MJD[1]; yb_shift$MJD[1]
common_mjd[length(common_mjd)]; sr_shift$MJD[length(sr_shift$MJD)]; yb_shift$MJD[length(yb_shift$MJD)]

min(sr_shift$MJD) < common_mjd[1]
min(yb_shift$MJD) < common_mjd[1]

max(sr_shift$MJD) > common_mjd[length(common_mjd)]
max(yb_shift$MJD) > common_mjd[length(common_mjd)]

################################################################################
#  Compute frequencies, create constants for correction
################################################################################

# frequency for Sr clock
compute_nuSr_ErYb <- function(data){
    nuSi = -105e6 + 388752 * (1e9 + data$SDR.frep_ErYb) / 2 - 100e6
    return(nuSr = (1716882 / 777577) * (nuSi - 216e6))
}

# frequency for Yb clock
compute_nuYb_ErYb <- function(data){
    nuYb = 2 * (-105e6 + 518237 * (1e9 + data$SDR.frep_ErYb) / 2 - data$SDR.fb_Yb_ErYb)
    return(nuYb)
}

nuSr <- compute_nuSr_ErYb(ErYb_data)
nuYb <- compute_nuYb_ErYb(ErYb_data)


YbSrRatio2020 <- 1.2075070393433378482
masercorrection <- -7.36631e-12
GR_shift_Yb <- -8.109e-16
GR_shift_Sr <- 10.660e-16
GR_shift_sea_level <- -1798.501e-16

total_correction_Yb <- 1 + GR_shift_Yb + GR_shift_sea_level + masercorrection
total_correction_Sr <- 1 + GR_shift_Sr + GR_shift_sea_level + masercorrection


################################################################################
## Analysis 1) Keep NaNs, concatenate, then linear interpolate
# Note: maybe can try reverse, interpolate then concatenate?
################################################################################

sr_shift1 <- sr_shift %>% filter(!is.na(shift))
yb_shift1 <- yb_shift %>% filter(!is.na(shift)) %>% select(-IS_GOOD)
#sr_shift1$shift <- format(sr_shift1$shift, digits=22)
#yb_shift1$shift <- format(yb_shift1$shift, digits=22)
#write.table(sr_shift1, file="sr_shift1_impute.csv", row.names=FALSE)
#write.table(yb_shift1, file="yb_shift1_impute.csv", row.names=FALSE)

min(sr_shift1$MJD) < common_mjd[390]
min(yb_shift1$MJD) < common_mjd[390]

max(sr_shift1$MJD) > common_mjd[length(common_mjd) - 8]
max(yb_shift1$MJD) > common_mjd[length(common_mjd) - 8]


shift_comb1 <- data.frame(sr_shift_comb = interp1(sr_shift1$MJD, sr_shift1$shift, xi=common_mjd[390:(length(common_mjd)-8)]),
                          yb_shift_comb = interp1(yb_shift1$MJD, yb_shift1$shift, xi=common_mjd[390:(length(common_mjd)-8)]))


nuSr_sub = nuSr[390:(length(common_mjd) - 8)]
nuYb_sub = nuYb[390:(length(common_mjd) - 8)]
shift_corrected1 <- data.frame(frequency_Sr_ErYb = (nuSr_sub + shift_comb1$sr_shift_comb) * total_correction_Sr,
                               frequency_Yb_ErYb = (nuYb_sub + shift_comb1$yb_shift_comb) * total_correction_Yb)

frequency_ratio_ErYb1 <- (shift_corrected1$frequency_Yb_ErYb/shift_corrected1$frequency_Sr_ErYb  - 
                            YbSrRatio2020) / YbSrRatio2020

plot.new()
p1 <- ggplot(data = data.frame(x=common_mjd[390:(length(common_mjd)-8)], y = frequency_ratio_ErYb1)) +
            geom_line(aes(x, y)) + geom_hline(yintercept=0, col="red", lwd=2)
p1

ggsave('freq_ratio1.tiff', p1, device = "tiff", width=10, height =5, dpi = 400) 
dev.off()

################################################################################
## Analysis 2) Two steps, linear imputation then linear interpolation
## other options include spline, stine
################################################################################
sr_shift2 <- na_interpolation(sr_shift$shift, option = "linear")
yb_shift2 <- na_interpolation(yb_shift$shift, option = "linear")
# sr_df2 <- data.frame(MJD = format(sr_shift$MJD, digits=15), shift = format(sr_shift2, digits=22))
# yb_df2 <- data.frame(MJD = format(yb_shift$MJD, digits=15), shift = format(yb_shift2, digits=22))
# write.table(sr_df2, file="sr_shift2_impute.csv", row.names=FALSE)
# write.table(yb_df2, file="yb_shift2_impute.csv", row.names=FALSE)

min(sr_shift$MJD) < common_mjd[1]
min(yb_shift$MJD) < common_mjd[1]

max(sr_shift$MJD) > common_mjd[length(common_mjd)]
max(yb_shift$MJD) > common_mjd[length(common_mjd)]

shift_comb2 <- data.frame(sr_shift_comb = interp1(sr_shift$MJD, sr_shift2, xi=common_mjd),
                          yb_shift_comb = interp1(yb_shift$MJD, yb_shift2, xi=common_mjd))

shift_corrected2 <-  data.frame(frequency_Sr_ErYb = (nuSr[1:22063] + shift_comb2$sr_shift_comb) * total_correction_Sr,
                               frequency_Yb_ErYb = (nuYb[1:22063] + shift_comb2$yb_shift_comb) * total_correction_Yb)

frequency_ratio_ErYb2 <- (shift_corrected2$frequency_Yb_ErYb/shift_corrected2$frequency_Sr_ErYb  - 
                            YbSrRatio2020) / YbSrRatio2020

plot.new()
p2 <- ggplot(data = data.frame(x=common_mjd, y = frequency_ratio_ErYb2)) +
            geom_line(aes(x, y)) + geom_hline(yintercept=0, col="red", lwd=2)
p2

ggsave('freq_ratio2.tiff', p2, device = "tiff", width=10, height =5, dpi = 400) 
dev.off()

################################################################################
## Analysis 3) Two steps, kalman filter imputation (imputeTS) then linear interpolation
## alternatives: na_locf (last obs carried forward), na_ma (moving average),
##               na_mean (mean imputation), na_random (random value with bound)
################################################################################
sr_shift3 <- na_kalman(sr_shift$shift, smooth = TRUE)
yb_shift3 <- na_kalman(yb_shift$shift, smooth = TRUE)
# sr_df3 <- data.frame(MJD = format(sr_shift$MJD, digits=15), shift = format(sr_shift3, digits=22))
# yb_df3 <- data.frame(MJD = format(yb_shift$MJD, digits=15), shift = format(yb_shift3, digits=22))
# write.table(sr_df3, file="sr_shift3_impute.csv", row.names=FALSE)
# write.table(yb_df3, file="yb_shift3_impute.csv", row.names=FALSE)

shift_comb3 <- data.frame(sr_shift_comb = interp1(sr_shift$MJD, sr_shift3, xi=common_mjd),
                          yb_shift_comb = interp1(yb_shift$MJD, yb_shift3, xi=common_mjd))

shift_corrected3 <-  data.frame(frequency_Sr_ErYb = (nuSr[1:22063] + shift_comb3$sr_shift_comb) * total_correction_Sr,
                               frequency_Yb_ErYb = (nuYb[1:22063] + shift_comb3$yb_shift_comb) * total_correction_Yb)

frequency_ratio_ErYb3 <- (shift_corrected3$frequency_Yb_ErYb/shift_corrected3$frequency_Sr_ErYb  - 
                            YbSrRatio2020) / YbSrRatio2020


plot.new()
p3 <- ggplot(data = data.frame(x=common_mjd, y = frequency_ratio_ErYb3)) +
            geom_line(aes(x, y)) + geom_hline(yintercept=0, col="red", lwd=2)
p3

ggsave('freq_ratio3.tiff', p3, device = "tiff", width=10, height =5, dpi = 400) 
dev.off()

################################################################################
################################################################################
## Scratch work
################################################################################
#For datetime manipulation: https://stackoverflow.com/questions/9839343/extracting-time-from-posixct
class(comb$MJD)
min(comb$MJD)
max(comb$MJD)

class(sr_shift$shift)
# #tst = sr_shift$shift
# tst[200]
# formatC(tst[200], digits=30, format="f")
# formatC(as.numeric(tst[200]), digits=30, format="f") %>% class
# formatC_tst = sapply(tst, formatC, digits=30, format="f")
# min(formatC_tst)

class(sr_shift$MJD)
min(sr_shift$MJD)
max(sr_shift$MJD)



# comb$time = convertToDateTime(comb$MJD, origin = "1858-11-17")
# summary(comb$MJD-60507)
# dtparts = t(as.data.frame(strsplit(comb$time,' ')))
# row.names(dtparts) = NULL
# 
# sr_shift$MJD = as.numeric(sr_shift$MJD, length=33) - 60507
# summary(sr_shift$MJD)
# sr_shift$time = convertToDateTime(sr_shift$MJD, origin = "1858-11-17")
# sr_shift$shift= as.numeric(sr_shift$shift, length=33)



ggplot_na_distribution(x=sr_shift$shift)#, x_axis_labels = sr_shift$MJD)
ggplot_na_distribution2(sr_shift$shift)
ggplot_na_gapsize(sr_shift$shift)

imp_sr = na_interpolation(sr_shift$shift)
ggplot_na_imputations(x_with_na=sr_shift$shift, x_with_imputations=imp_sr) #, x_axis_labels = sr_shift$MJD)

imp_yb = na_interpolation(yb_shift$shift)
ggplot_na_imputations(x_with_na=yb_shift$shift, x_with_imputations=imp_yb) #, x_axis_labels = yb_shift$MJD)


# ## Add column of MJD differences
# MJD_diff = rep(NA, length(comb$MJD))
# for(i in 1:(length(comb$MJD)-1)){
#   MJD_diff[i+1] = as.double(comb$MJD[i+1])-as.double(comb$MJD[i])
# }
# print(summary(MJD_diff), digits=22)
# comb <- comb %>% mutate(diff = MJD_diff)

# MJD_diff = rep(NA, length(sr_shift$MJD))
# for(i in 1:(length(sr_shift$MJD)-1)){
#   MJD_diff[i+1] = as.double(sr_shift$MJD[i+1])-as.double(sr_shift$MJD[i])
# }
# sr_shift <- sr_shift %>% mutate(diff = MJD_diff)

# MJD_diff = rep(NA, length(yb_shift$MJD))
# for(i in 1:(length(yb_shift$MJD)-1)){
#   MJD_diff[i+1] = as.double(yb_shift$MJD[i+1])-as.double(yb_shift$MJD[i])
# }
# yb_shift <- yb_shift %>% mutate(diff = MJD_diff)

## Matching MJD values
sr_shift$MJD[1]
yb_shift$MJD[1]
comb$MJD[1]; sr_shift$MJD[325]; yb_shift[140]
comb$MJD[22064]; yb_shift$MJD[25396]
comb$MJD[22113]; sr_shift$MJD[3132]
comb$MJD[24289]

