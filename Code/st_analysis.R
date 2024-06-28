library("openxlsx")
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

ErYb_AlYb = read.csv2("Data/clockRatioTimeseries_first2024dryRun/2024-04-24-ErYb-AlYb.csv",
                      header = FALSE,sep="",na.strings = "NaN")
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

ggplot_na_distribution(x=ErYb_AlYb$f_hz_num_shift, x_axis_labels = ErYb_AlSr$date) 
summary(ErYb_AlYb$f_hz_num_shift)

ggplot_na_distribution(x=ErYb_YbSr$f_hz_num_shift, x_axis_labels = ErYb_AlSr$date) 
summary(ErYb_YbSr$f_hz_num_shift)
