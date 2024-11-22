# -*- coding: utf-8 -*-
import pandas as pd
from decimal import Decimal
import matplotlib.pyplot as plt
import numpy as np
 
 
################################################################################
#############  Functions for data loading #####################################
################################################################################

## comb data 
def open_ErYb_data(data_path):
    # keys to read out as string
    key2read = ["fo_ErYb", "SDR:fb_Yb_ErYb", "SDR:frep_ErYb", "fb_Si_ErYb"]
    types = {key: str for key in key2read}
    types["MJD"] = float
 
    # # Read the CSV file
    data = pd.read_csv(data_path, header=1, delimiter="\t", dtype=types, engine="python")
 
    # Convert the strings to Decimal for the given keys
    #Q: Why not do the same when importing shift data? 
    #A: not necessary to do here
    for k in key2read:
        data[k] = data[k].apply(Decimal)
 
    # reindex data
    data.index = range(len(data))
 
    return data[list(types.keys())]
 
## Sr shift data 
def open_shiftfile_Sr(datapath):
    #Q: What does dtype={1: str} mean? Nick unsure... look into. string to float to decimal. 
    data = pd.read_csv(datapath, header=21, delimiter="\t", dtype={1: str}, engine="python")
 
    # Replace column names
    data.columns = ["MJD", "shift", "IS_GOOD"]
 
    # Change column type from float to bool
    #Note: lambda is way to define small local function definition 
    data["IS_GOOD"] = data["IS_GOOD"].apply(lambda x: x == 1.0)
 
    # Put NaN in data["shift"] where data["IS_GOOD"] is 0
    data.loc[~data["IS_GOOD"], "shift"] = np.nan
 
    # Change column type to float
    data["shift"] = data["shift"].apply(float)
 
    return data
 
 ## Yb shift data
def open_shiftfile_Yb(datapath):
    data = pd.read_csv(datapath, header=7, delimiter=r"\t", dtype={1: str}, engine="python")
 
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
 
# frequency for Sr clock
#Note: dependent upon silicon cavity, not a clock just a laser 
def compute_nuSr_ErYb(data):
    data["nuSi"] = -Decimal("105e6") + Decimal("388752") * (Decimal("1e9") + data["SDR:frep_ErYb"]) / Decimal(2) - Decimal("100e6")
    data["nuSr"] = (Decimal("1716882") / Decimal("777577")) * (data["nuSi"] - Decimal("216e6"))


# freuency for Yb clock
def compute_nuYb_ErYb(data):
    data["nuYb"] = -Decimal("105e6") + Decimal("518237") * (Decimal("1e9") + data["SDR:frep_ErYb"]) / Decimal(2) - data["SDR:fb_Yb_ErYb"]
    data["nuYb"] = Decimal(2) * data["nuYb"]
 
 
 
 
 
################################################################################
#############################  Load data #######################################
################################################################################
 
# this is the path where all data is stored on my computer
#path = "/Users/nvn/Box Sync/NIST/data/2024 clock comparison data/20240716/"
path = "/Users/smt3/Documents/test-python/"

# load comb data
data_ErYb = open_ErYb_data(path + "20240716_Deglitched_ErYb_only.dat")
 
# load Sr shift data
shift_data_Sr = pd.concat(
    [
        open_shiftfile_Sr(path + "20240716_Sr_Freq_Shifts_1.dat"),
        open_shiftfile_Sr(path + "20240716_Sr_Freq_Shifts_2.dat"),
        open_shiftfile_Sr(path + "20240716_Sr_Freq_Shifts_3.dat"),
        open_shiftfile_Sr(path + "20240716_Sr_Freq_Shifts_4.dat"),
    ],
    ignore_index=True,
)
 
# load Yb shift data
shift_data_Yb = open_shiftfile_Yb(path + "20240716_Yb_Freq_Shifts.txt")
 
 
 
 
 
################################################################################
###############  get Sr and Yb optical frequencies #############################
################################################################################
 
compute_nuSr_ErYb(data_ErYb)
compute_nuYb_ErYb(data_ErYb)
 
 
################################################################################
#########################  Interpolation #######################################
################################################################################
 
print("Sr sampling rate", '{:0.2}'.format( np.mean(np.diff(shift_data_Sr["MJD"]))*24*3600 ), "s" )
print("Yb sampling rate", '{:0.2}'.format( np.mean(np.diff(shift_data_Yb["MJD"]))*24*3600 ), "s" )
 
common_mjd = data_ErYb["MJD"]
nuYb = data_ErYb["nuYb"]
nuSr = data_ErYb["nuSr"]
 
 
# look for gaps in Sr shift file longer than 15s and fill them with NaN
interval1s = 15/(3600*24)
MJDs1s = []   
shifts1s = []
for i in range(0, len(shift_data_Sr["MJD"])):
    newMJD = shift_data_Sr["MJD"][i]
    MJDs1s.append(newMJD)
    shifts1s.append(shift_data_Sr["shift"][i])
    
    if (i+1) < len(shift_data_Sr["MJD"]):
        while abs(newMJD - shift_data_Sr["MJD"][i+1]) > interval1s:
            newMJD = newMJD + interval1s
            MJDs1s.append(newMJD)
            shifts1s.append(np.nan)
 
 
del shift_data_Sr
shift_data_Sr = pd.DataFrame({})
shift_data_Sr["MJD"] = MJDs1s
shift_data_Sr["shift"] = shifts1s
 
 
# interpolate
# np package function for interpolation: https://numpy.org/doc/stable/reference/generated/numpy.interp.html 
# np.interp(<want to obs x>, <observed x>, <observed y>)
shiftSr = np.interp(data_ErYb["MJD"], shift_data_Sr["MJD"], shift_data_Sr["shift"])
shiftYb = np.interp(data_ErYb["MJD"], shift_data_Yb["MJD"], shift_data_Yb["shift"])
 
# change type to high-precision Decimal
# Q: why isn't this done at the very beginning, when data is read in? Or at least before interpolation?
# A: b/c np.interp works on floats only. is faster this way. 
shiftSr = [Decimal(i) for i in shiftSr]
shiftYb = [Decimal(i) for i in shiftYb]
 
 
 
################################################################################
#############################  Plotting  #######################################
################################################################################
 
# Yb/Sr ratio from 2020 BACON paper
YbSrRatio2020 = Decimal("1.2075070393433378482") ##underlying "truth"
 
# frequency corrections
masercorrection = Decimal("-7.36631e-12")
GR_shift_Yb = Decimal("-8.109e-16")
GR_shift_Sr = Decimal("10.660e-16")
GR_shift_sea_level = Decimal("-1798.501e-16")
 
total_correction_Yb = Decimal("1") + GR_shift_Yb + GR_shift_sea_level + masercorrection
total_correction_Sr = Decimal("1") + GR_shift_Sr + GR_shift_sea_level + masercorrection
 
# add comb frequencies and clock shift files
frequency_Yb_ErYb = [(i + j) * total_correction_Yb for i,j in zip(nuYb, shiftYb)]
frequency_Sr_ErYb = [(i + j) * total_correction_Sr for i,j in zip(nuSr, shiftSr)]
 
# Yb/Sr ratio offset from BACON paper
#Q: Why dividing the offset by the YbSrRatio2020? 
#A: convention, normalize to one. ideally answer would be 1
frequency_ratio_ErYb = [(i / j - YbSrRatio2020)/YbSrRatio2020 for i,j in zip(frequency_Yb_ErYb, frequency_Sr_ErYb)]
#Note: Is nanmean necessary? after interpolation? still missing data after interpolation but what does missingess look like before/after interpolation
# Nick note: do not want to interpolate across the long periods of missing data, missing data "should be" from experimental issues 
# Note contd: goal is for shift files to be grid of points 1 sec apart. look for weird patterns in missingness  
print("Yb/Sr ratio offset from BACON paper", '{:0.5}'.format(np.nanmean(frequency_ratio_ErYb)) )
 
plt.figure()
plt.plot(common_mjd, frequency_ratio_ErYb, '.')
plt.grid(visible=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.xlabel("MJD")
plt.ylabel("offset from Yb/Sr ratio BACON ratio")
plt.title("Yb/Sr ratio from July 16 2024")
plt.show()
