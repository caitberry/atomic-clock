library("openxlsx")
library(tidyverse)
##-------------------------------------------------------------------------------------------
## Shift data files 
##-------------------------------------------------------------------------------------------
Al_shift = read.delim("Data/20240426_Al_Freq_Shifts_ErYb.dat",header = TRUE,
                        skip=21, sep="") #,na.strings = "NaN")
Yb_shift = read.delim("Data/20240426_Yb_Freq_Shifts.txt",header = TRUE,
                        skip=7, sep="")#,na.strings = "NaN")
Sr_shift = read.delim("Data/20240426_clock_lock.dat",header = TRUE,
                        skip=21, sep="")#,na.strings = "NaN")


Al_shift = Al_shift %>% mutate(missing=is.na(f_shift.Hz.),
                                date=convertToDateTime(t_mjd, origin = "1858-11-17"))
Yb_shift = Yb_shift %>% mutate(missing=is.na(f_shift.Hz.),
                               date=convertToDateTime(t_mjd, origin = "1858-11-17"))
Sr_shift = Sr_shift %>% mutate(missing=is.na(f_shift.Hz.),
                               date=convertToDateTime(t_mjd, origin = "1858-11-17"))

## No missing data in these files 
Yb_shift$missing %>% sum

##-------------------------------------------------------------------------------------------
## clockRatioTimeseries_first2024dryRun data 
##-------------------------------------------------------------------------------------------
## Date:4/24/24 
ErYb_AlSr = read.csv2("Data/clockRatioTimeseries_first2024dryRun/2024-04-24-ErYb-AlSr.csv",
           header = FALSE,sep="",na.strings = "NaN")
names(ErYb_AlSr) = c("t_mjd", "f_hz")
ErYb_AlSr = ErYb_AlSr %>% mutate(missing = is.na(f_hz),
                                 date = convertToDateTime(t_mjd, origin = "1858-11-17"),
                                 f_hz_num = as.numeric(f_hz))  %>% ##22 digits 
                          mutate(f_hz_num_shift = f_hz_num - min(na.omit(f_hz_num)))
ErYb_AlSr$missing %>% sum

#ErYb_AlYb = read.csv2("Data/clockRatioTimeseries_first2024dryRun/2024-04-24-ErYb-AlYb.csv",
#                      header = FALSE,sep="",na.strings = "NaN")
ErYb_AlYb = read.csv2("Data/clockRatioOffsetTimeseries_first2024dryRun/2024-04-24-ErYb-AlYb.csv",
                      header = FALSE,skip=3,sep="",na.strings = "nan")
names(ErYb_AlYb) = c("t_mjd", "f_hz")
ErYb_AlYb = ErYb_AlYb %>% mutate(missing = is.na(f_hz),
                                 date = convertToDateTime(t_mjd, origin = "1858-11-17"),
                                 f_hz_num = as.numeric(f_hz))  %>% ##22 digits 
                          mutate(f_hz_num_shift = f_hz_num - min(na.omit(f_hz_num)))
ErYb_AlYb$missing %>% sum

ErYb_YbSr = read.csv2("Data/clockRatioTimeseries_first2024dryRun/2024-04-24-ErYb-YbSr.csv",
                      header = FALSE,sep="",na.strings = "NaN")
names(ErYb_YbSr) = c("t_mjd", "f_hz")
ErYb_YbSr = ErYb_YbSr %>% mutate(missing = is.na(f_hz),
                                 date = convertToDateTime(t_mjd, origin = "1858-11-17"),
                                 f_hz_num = as.numeric(f_hz))  %>% ##22 digits 
                          mutate(f_hz_num_shift = f_hz_num - min(na.omit(f_hz_num)))
ErYb_YbSr$missing %>% sum

## Visualize the time series with missing data
## Reference: https://cran.r-project.org/web/packages/imputeTS/vignettes/gallery_visualizations.html
library("tidyverse")
library("imputeTS")
ggplot_na_distribution(x=ErYb_AlSr$f_hz_num_shift, x_axis_labels = ErYb_AlSr$date) 
summary(ErYb_AlSr$f_hz_num_shift)
ggplot_na_distribution2(ErYb_AlSr$f_hz_num_shift)
ggplot_na_gapsize(ErYb_AlSr$f_hz_num_shift)
imp = na_interpolation(ErYb_AlSr$f_hz_num_shift)
ggplot_na_imputations(x_with_na=ErYb_AlSr$f_hz_num_shift, x_with_imputations=imp, x_axis_labels = ErYb_AlSr$date)

tmp = ErYb_AlYb[!ErYb_AlYb$missing,]
ggplot_na_distribution(x=tmp$f_hz_num_shift, x_axis_labels = tmp$date) 
ggplot_na_distribution(x=ErYb_AlYb$f_hz_num_shift, x_axis_labels = ErYb_AlSr$date) 
summary(ErYb_AlYb$f_hz_num_shift)

ggplot_na_distribution(x=ErYb_YbSr$f_hz_num_shift, x_axis_labels = ErYb_AlSr$date) 
summary(ErYb_YbSr$f_hz_num_shift)

##--------------------------------------------------------------------------
## Note: Offset and original data are not the same length 
##       also - offset file has converted f_hz into a character representing a number in scientific notation (e.g. 5.343434e-15)
##---------------------------------------------------------------------------
library(tidyverse); library(readr)
#options(digits=22)
#readr::parse_double(), sprintf(), as.numeric(), formatC(), format()
ErYb_AlYb_offset = read.csv2("Data/clockRatioOffsetTimeseries_first2024dryRun/2024-04-24-ErYb-AlYb.csv",
                      header = FALSE, numerals = "no.loss", #colClasses=c("character", "numeric"), 
                      skip=3,sep="",na.strings = "nan")
ErYb_AlYb = read.csv2("Data/clockRatioTimeseries_first2024dryRun/2024-04-24-ErYb-AlYb.csv",
                             header = FALSE, numerals="no.loss", #colClasses=c("character", "numeric"), 
                      sep="",na.strings = "NaN")
full_offset = ErYb_AlYb_offset$V2[!is.na(ErYb_AlYb_offset$V2)]
full_orig = ErYb_AlYb$V2[!is.na(ErYb_AlYb$V2)]
length(full_offset)
length(full_orig)
plot(full_offset)
plot(full_orig)
full_offset[100]
full_orig[100]
sprintf("%s", full_orig[100]) %>% class()

(formatC(full_orig[100],digits=30)-as.character(2.162887127516663703))/2.162887127516663703
##figure out if one is a subset of the other, try adding NAs at beginning in end or comparing size of non missing data

ErYb_AlYb = ErYb_AlYb %>% #mutate(t_mjd = parse_double(V1), f_hz=sprintf(V2, fmt='%.30f')) %>% 
                          mutate(date = convertToDateTime(V1, origin = "1858-11-17"))
ggplot_na_distribution(x=ErYb_AlYb$V2, x_axis_labels = ErYb_AlYb$date)
plot(y=ErYb_AlYb$V2, x = ErYb_AlYb$date)

## Q: how to format the data to analyze



