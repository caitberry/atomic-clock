import pandas as pd
from decimal import Decimal
import matplotlib.pyplot as plt
import numpy as np
 
 
################################################################################
#############  Functions for data loading #####################################
################################################################################

## comb data 
def open_ErYb_data(data_path, header=2):
    # keys to read out as string
    key2read = ["MJD", "timer", "SDR:frep_ErYb", "fo_ErYb", "fb_Si_ErYb", "fb_Al_ErYb", "fb_Yb_ErYb"] 
    types = {key: str for key in key2read}
    types["MJD"] = float
 
    # # Read the CSV file
    data = pd.read_csv(data_path, header=1, delimiter="\t", dtype=types, engine="python")
 
    # Convert the strings to Decimal for the given keys
    for k in key2read:
        data[k] = data[k].apply(Decimal)
 
    # reindex data
    data.index = range(len(data))
 
    return data[list(types.keys())]

## Al shift data 
def open_shiftfile_Al(datapath):
    data = pd.read_csv(datapath, header=30, delimiter="\t", dtype={1: str}, engine="python")
 
    # Replace column names
    data.columns = ["MJD", "shift", "IS_GOOD"]
 
    # Change column type from float to bool
    data["IS_GOOD"] = data["IS_GOOD"].apply(lambda x: x == 1.0)
 
    # Put NaN in data["shift"] where data["IS_GOOD"] is 0
    data.loc[~data["IS_GOOD"], "shift"] = np.nan
 
    # Change column type to float
    data["shift"] = data["shift"].apply(float)
 
    return data
 
## Sr shift data 
def open_shiftfile_Sr(datapath):
    data = pd.read_csv(datapath, header=22, delimiter="\t", dtype={1: str}, engine="python")
 
    # Replace column names
    data.columns = ["MJD", "shift", "IS_GOOD"]
 
    # Change column type from float to bool
    data["IS_GOOD"] = data["IS_GOOD"].apply(lambda x: x == 1.0)
 
    # Put NaN in data["shift"] where data["IS_GOOD"] is 0
    data.loc[~data["IS_GOOD"], "shift"] = np.nan
 
    # Change column type to float
    data["shift"] = data["shift"].apply(float)
 
    return data
 
 ## Yb shift data
def open_shiftfile_Yb(datapath):
    data = pd.read_csv(datapath, header=8, delimiter=r"\t", dtype={1: str}, engine="python")
 
    # Replace column names
    data.columns = ["MJD", "shift", "IS_GOOD"]
 
    # Change column type from float to bool
    data["IS_GOOD"] = data["IS_GOOD"].apply(lambda x: x == 1.0)
 
    # Put NaN in data["shift"] where data["IS_GOOD"] is 0
    data.loc[~data["IS_GOOD"], "shift"] = np.nan
 
    # Change column type to float
    data["shift"] = data["shift"].apply(float)
 
    return data
 
 
 
 
################################################################################
#############  Functions to find optical frequencies with comb equation ########
################################################################################
 
# frequency for Al+ clock
def compute_nuAl_ErYb(data):
    data["nuAl"] = -Decimal("105e6") + Decimal("560444") * (Decimal("1e9") + data["SDR:frep_ErYb"]) / Decimal(2) - data["fb_Al_ErYb"]
    data["nuAl"] = Decimal(4) * data["nuAl"]   

# frequency for Sr clock 
def compute_nuSr_ErYb(data):
    data["nuSi"] = -Decimal("105e6") + Decimal("388752") * (Decimal("1e9") + data["SDR:frep_ErYb"]) / Decimal(2) - Decimal("100e6")
    data["nuSr"] = (Decimal("1716882") / Decimal("777577")) * (data["nuSi"] - Decimal("216e6"))


# freuency for Yb clock
def compute_nuYb_ErYb(data):
    data["nuYb"] = -Decimal("105e6") + Decimal("518237") * (Decimal("1e9") + data["SDR:frep_ErYb"]) / Decimal(2) - data["fb_Yb_ErYb"]
    data["nuYb"] = Decimal(2) * data["nuYb"]
 
 
 
################################################################################
#############################  Load data #######################################
################################################################################
 
path = "/Users/smt3/Documents/GitHub/atomic-clock/st-interp/three_clocks/"

# load comb data
data_ErYb = open_ErYb_data(path + "20240813_Deglitched_ErYb_only1.dat")
 
# load Al shift data 
shift_data_Al = open_shiftfile_Al(path + "20240813_Al+_Freq_Shifts_ErYb.dat")

# load Sr shift data
shift_data_Sr = open_shiftfile_Sr(path + "20240813_Sr_Freq_Shifts.dat")
 
# load Yb shift data
shift_data_Yb = open_shiftfile_Yb(path + "20240813_Yb_Freq_Shifts.txt")
 
 
 
 
 
################################################################################
###############  get Sr and Yb optical frequencies #############################
################################################################################
 
compute_nuSr_ErYb(data_ErYb)
compute_nuYb_ErYb(data_ErYb)
compute_nuAl_ErYb(data_ErYb)

################################################################################
## Visualize gaps in data 
################################################################################
## pip install missingno 

import missingno as msno
atts = ["MJD", "shift"]

miss_inx_Al = pd.isnull(shift_data_Al["shift"])
print(miss_inx_Al.sum())
msno.matrix(shift_data_Al[atts])
plt.show()

miss_inx_Sr = pd.isnull(shift_data_Sr["shift"])
print(miss_inx_Sr.sum())
msno.matrix(shift_data_Sr[atts])
plt.show()

miss_inx_Yb = pd.isnull(shift_data_Yb["shift"])
print(miss_inx_Yb.sum())
msno.matrix(shift_data_Yb[atts])
plt.show()
 
 
################################################################################
#########################  Data Processing #####################################
################################################################################

## Extract only "IS_GOOD" data for analysis------------------
good_condition_al = shift_data_Al["IS_GOOD"] == 1
shift_data_Al_good = shift_data_Al[good_condition_al].reset_index(drop=True)
good_condition_sr = shift_data_Sr["IS_GOOD"] == 1
shift_data_Sr_good = shift_data_Sr[good_condition_sr].reset_index(drop=True)
good_condition_yb = shift_data_Yb["IS_GOOD"] == 1
shift_data_Yb_good = shift_data_Yb[good_condition_yb].reset_index(drop=True)


## Find common MDJ values--------------------------------
#NOTE: the following does not take into account adjustements made based on large chunks of missing data

#Change comb mjd type to float  
common_mjd = data_ErYb["MJD"].astype(float)
nuAl = data_ErYb["nuAl"].astype(float)
nuYb = data_ErYb["nuYb"].astype(float)
nuSr = data_ErYb["nuSr"].astype(float)
#convert comb nu values to high precision decimal 
nuAl = [Decimal(i) for i in nuAl]
nuSr = [Decimal(i) for i in nuSr]
nuYb = [Decimal(i) for i in nuYb] 

# Length of the 'MJD' column
len_comb = len(common_mjd) 
len_Al = len(shift_data_Al_good['MJD'])                  
len_Sr = len(shift_data_Sr_good['MJD'])
len_Yb = len(shift_data_Yb_good['MJD'])

#function to extract element as close to target as possible w/out going over
def lb_extract(target, data):
    inx = 0
    stopper = 1
    while stopper == 1:
        if data[inx] <= target:
            inx += 1
        else:
            return inx  

#function to extract element as close to target as possible w/out going under 
def ub_extract(target, data):
    inx = 1
    stopper = 1
    while stopper == 1:
        if data[len(data)-inx] >= target:
            inx += 1
        else:
            return len(data)-inx  

#Compare start and end points (assuming no missing data b/c IS_GOOD variable already accounted for)
print("first comb time point: ", common_mjd[0])
print("first good Al time point: ", shift_data_Al_good["MJD"][0])
print("first good Sr time point: ", shift_data_Sr_good["MJD"][0])
print("first good Yb time point: ", shift_data_Yb_good["MJD"][0])
start_times = {common_mjd[0], shift_data_Al_good["MJD"][0], shift_data_Sr_good["MJD"][0], shift_data_Yb_good["MJD"][0]}
print("----->Latest start time: ", max(start_times)) 
print("last comb time point: ", common_mjd[len_comb-1])
print("last good Al time point: ", shift_data_Al_good["MJD"][len_Al-1])
print("last good Sr time point: ", shift_data_Sr_good["MJD"][len_Sr-1])
print("last good Yb time point: ", shift_data_Yb_good["MJD"][len_Yb-1])
end_times = {common_mjd[len_comb-1], shift_data_Al_good["MJD"][len_Al-1], shift_data_Sr_good["MJD"][len_Sr-1], shift_data_Yb_good["MJD"][len_Yb-1]}
print("----->Earliest end time: ", min(end_times), "\n")


#NOTE: The following must be edited based on the resuts of the print statements above
print("since Sr starts the latest, start comb observations at index ", lb_extract(target = shift_data_Sr_good["MJD"][0], data = common_mjd), 
         ", start Al observations at index ", lb_extract(target = shift_data_Sr_good["MJD"][0], data = shift_data_Al_good["MJD"]),
         " and start Yb observations at index ", lb_extract(target = shift_data_Sr_good["MJD"][0], data = shift_data_Yb_good["MJD"]))
print("since Yb ends first, end comb observations at index ", ub_extract(target = shift_data_Yb_good["MJD"][len_Yb-1], data = common_mjd), 
         ", end Sr observations at index ", ub_extract(target = shift_data_Yb_good["MJD"][len_Yb-1], data = shift_data_Sr_good["MJD"]),
         " and end Al observations at index ", ub_extract(target = shift_data_Yb_good["MJD"][len_Yb-1], data = shift_data_Al_good["MJD"]), "\n")

#comb MJD index 
comb = pd.DataFrame()
comb_init = lb_extract(target = shift_data_Sr_good["MJD"][0], data = common_mjd) #edit
comb_end = ub_extract(target = shift_data_Yb_good["MJD"][len_Yb-1], data = common_mjd) #edit
comb["MJD"] = common_mjd[comb_init:comb_end]
comb["nuAl"] = nuAl[comb_init:comb_end]
comb["nuYb"] = nuYb[comb_init:comb_end]
comb["nuSr"] = nuSr[comb_init:comb_end]

#Al MJD index 
shift_data_Al = shift_data_Al_good[lb_extract(target = shift_data_Sr_good["MJD"][0], data = shift_data_Al_good["MJD"]):ub_extract(target = shift_data_Yb_good["MJD"][len_Yb-1], data = shift_data_Al_good["MJD"])]
#Sr MJD index 
shift_data_Sr = shift_data_Sr_good[0:ub_extract(target = shift_data_Yb_good["MJD"][len_Yb-1], data = shift_data_Sr_good["MJD"])]
#Yb MJD index  
shift_data_Yb = shift_data_Yb_good[lb_extract(target = shift_data_Sr_good["MJD"][0], data = shift_data_Yb_good["MJD"]):len_Yb-1]


################################################################################
#########################  Interpolation techniques  ###########################
################################################################################

shiftAl = np.interp(common_mjd, shift_data_Al["MJD"], shift_data_Al["shift"])
shiftSr = np.interp(common_mjd, shift_data_Sr["MJD"], shift_data_Sr["shift"])
shiftYb = np.interp(common_mjd, shift_data_Yb["MJD"], shift_data_Yb["shift"])

shiftAl = [Decimal(i) for i in shiftAl]
shiftSr = [Decimal(i) for i in shiftSr]
shiftYb = [Decimal(i) for i in shiftYb] 
 
################################################################################
#############################  Plotting  #######################################
################################################################################
 
# Ratios from 2020
YbSrRatio2020 = Decimal("1.2075070393433378482") 
AlYbRatio2020 = Decimal("2.162887127516663703")
AlSrRatio2020 = Decimal("2.611701431781463025")
 
# frequency corrections
masercorrection = Decimal("-7.36631e-12")
GR_shift_Al = Decimal("-8.114e-16")
GR_shift_Yb = Decimal("-8.109e-16")
GR_shift_Sr = Decimal("10.660e-16")
GR_shift_sea_level = Decimal("-1798.501e-16")

total_correction_Yb = Decimal("1") + GR_shift_Yb + GR_shift_sea_level + masercorrection
total_correction_Sr = Decimal("1") + GR_shift_Sr + GR_shift_sea_level + masercorrection
total_correction_Al = Decimal("1") + GR_shift_Al + GR_shift_sea_level + masercorrection
 
# add comb frequencies and clock shift files
frequency_Yb_ErYb = [(i + j) * total_correction_Yb for i,j in zip(nuYb, shiftYb)]
frequency_Sr_ErYb = [(i + j) * total_correction_Sr for i,j in zip(nuSr, shiftSr)]
frequency_Al_ErYb = [(i + j) * total_correction_Al for i,j in zip(nuAl, shiftAl)]
 
# Yb/Sr ratio offset from BACON paper
frequency_ratio_ErYb1 = [(i / j - YbSrRatio2020)/YbSrRatio2020 for i,j in zip(frequency_Yb_ErYb, frequency_Sr_ErYb)]
print("Yb/Sr ratio offset from BACON paper", '{:0.5}'.format(np.nanmean(frequency_ratio_ErYb1)) )
plt.figure()
plt.plot(common_mjd, frequency_ratio_ErYb1, '.')
plt.grid(visible=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.xlabel("MJD")
plt.ylabel("offset from Yb/Sr ratio BACON ratio")
plt.title("Yb/Sr ratio from 8/13/24 on common MJD grid")
plt.show()


# Al/Yb ratio offset 
frequency_ratio_ErYb2 = [(i / j - AlYbRatio2020)/AlYbRatio2020 for i,j in zip(frequency_Al_ErYb, frequency_Yb_ErYb)]
print("Al/Yb ratio offset from BACON paper", '{:0.5}'.format(np.nanmean(frequency_ratio_ErYb2)) )
plt.figure()
plt.plot(common_mjd, frequency_ratio_ErYb2, '.')
plt.grid(visible=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.xlabel("MJD")
plt.ylabel("offset from Al/Yb ratio BACON ratio")
plt.title("Al/Yb ratio from 8/13/24 on common MJD grid")
plt.show()


# Al/Sr ratio offset  
frequency_ratio_ErYb3 = [(i / j - AlSrRatio2020)/AlSrRatio2020 for i,j in zip(frequency_Al_ErYb, frequency_Sr_ErYb)]
print("Al/Sr ratio offset from BACON paper", '{:0.5}'.format(np.nanmean(frequency_ratio_ErYb3)) )
plt.figure()
plt.plot(common_mjd, frequency_ratio_ErYb3, '.')
plt.grid(visible=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.xlabel("MJD")
plt.ylabel("offset from Al/Sr ratio BACON ratio")
plt.title("Al/Sr ratio from 8/13/24 on common MJD grid")
plt.show()