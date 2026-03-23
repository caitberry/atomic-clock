#%% 
import pandas as pd
from decimal import Decimal
import matplotlib.pyplot as plt
import numpy as np
import math 
from typing import List, Union, Optional
from astropy.time import Time
import missingno as msno
import statsmodels.api as sm
import allantools 
from tabulate import tabulate ##new

def open_ErYb_data(data_path):
    key2read = ["MJD", "timer", "SDR:frep_ErYb", "fo_ErYb", "fb_Si_ErYb", "fb_Yb_ErYb", "fb_Al_ErYb"] 
    types = {key: str for key in key2read}
    types["MJD"] = float
 
    data = pd.read_csv(data_path, header=1, delimiter="\t", dtype=types, engine="python")
 
    for k in key2read:
        data[k] = data[k].apply(Decimal)
 
    data.index = range(len(data))
 
    return data[list(types.keys())]

def open_shiftfile_Al(datapath):
    data = pd.read_csv(datapath, header=30, delimiter="\t", dtype={1: str}, engine="python")
 
    data.columns = ["MJD", "shift", "IS_GOOD"]
 
    data["IS_GOOD"] = data["IS_GOOD"].apply(lambda x: x == 1.0)
 
    data.loc[~data["IS_GOOD"], "shift"] = np.nan
 
    data["shift"] = data["shift"].apply(float)
 
    return data

def open_shiftfile_Sr(datapath): ##Note: assuming (unverified) that for days_irregular all Sr .dat files have this format.... 
    data = pd.read_csv(datapath, header=24, delimiter="\t", dtype={1: str}, engine="python")

    data.columns = ["MJD", "shift", "IS_GOOD"]
 
    data["IS_GOOD"] = data["IS_GOOD"].apply(lambda x: x == 1.0)
 
    data.loc[~data["IS_GOOD"], "shift"] = np.nan
 
    data["shift"] = data["shift"].apply(float)
 
    return data
 
def open_shiftfile_Yb(datapath):
    data = pd.read_csv(datapath, header=8, delimiter=r"\t",  dtype={1: str}, engine="python")
 
    data.columns = ["MJD", "shift", "IS_GOOD"]
 
    data["IS_GOOD"] = data["IS_GOOD"].apply(lambda x: x == 1.0)
 
    data.loc[~data["IS_GOOD"], "shift"] = np.nan
 
    data["shift"] = data["shift"].apply(float)
 
    return data
 
def open_maser_correction(datapath):
    data = pd.read_csv(datapath, header=1, delimiter=",", dtype={1: str}, engine="python")
 
    data.columns = ["date", "maser_offset"]
 
    data["date"] = data["date"].apply(str)#.str.split("-").str.join("")
    data["maser_offset"] = data["maser_offset"].apply(float)
 
    return data

## ? YbSr ratio for all 13 days 
## ? AlYb and AlSr ratio only for 9 of the 13 days 
path = "/Users/smt3/Documents/GitHub/2025 clock comparison data/"
days = [20250116, 20250124, 20250204, 20250227, 20250304, 20250307, 20250318]
days = list(map(str, days))
day_index = 4 ##ignore 0 and 1 b/c is 2024 data; 2:6

days_irregular = [20250206, 20250228, 20250306, 20250313, 20250320, 20250321]
days_irregular = list(map(str, days_irregular))
day_irregular_index = None ##set to None if analyzing days, set to 0:5 if analyzing days_irregular

## For days_irregular, read in Sr data with the following 
if day_irregular_index != None:
    if days_irregular[day_irregular_index] == "20250206":
        shift_data_Sr = pd.concat(
            [
                open_shiftfile_Sr(path + days_irregular[day_irregular_index] + "/" + days_irregular[day_irregular_index] + "_clock_lock0.dat"),
                open_shiftfile_Sr(path + days_irregular[day_irregular_index] + "/" + days_irregular[day_irregular_index] + "_clock_lock2.dat"),
            ],
            ignore_index=True,
        )
    elif days_irregular[day_irregular_index] == "20250228":
        shift_data_Sr = open_shiftfile_Sr(path + days_irregular[day_irregular_index] + "/" + days_irregular[day_irregular_index] + "_clock_lock1.dat")
    elif days_irregular[day_irregular_index] == "20250306":
        shift_data_Sr = open_shiftfile_Sr(path + days_irregular[day_irregular_index] + "/" + days_irregular[day_irregular_index] + "_clock_lock2.dat")
    elif days_irregular[day_irregular_index] == "20250313":
        shift_data_Sr = pd.concat(
            [
                open_shiftfile_Sr(path + days_irregular[day_irregular_index] + "/" + days_irregular[day_irregular_index] + "_clock_lock0.dat"),
                open_shiftfile_Sr(path + days_irregular[day_irregular_index] + "/" + days_irregular[day_irregular_index] + "_clock_lock1.dat"),
                open_shiftfile_Sr(path + days_irregular[day_irregular_index] + "/" + days_irregular[day_irregular_index] + "_clock_lock2.dat"),
                open_shiftfile_Sr(path + days_irregular[day_irregular_index] + "/" + days_irregular[day_irregular_index] + "_clock_lock3.dat"),
            ],
            ignore_index=True,
        )
    elif days_irregular[day_irregular_index] == "20250320":
        shift_data_Sr = pd.concat(
            [
                open_shiftfile_Sr(path + days_irregular[day_irregular_index] + "/" + days_irregular[day_irregular_index] + "_clock_lock0.dat"),
                open_shiftfile_Sr(path + days_irregular[day_irregular_index] + "/" + days_irregular[day_irregular_index] + "_clock_lock2.dat"),
                open_shiftfile_Sr(path + days_irregular[day_irregular_index] + "/" + days_irregular[day_irregular_index] + "_clock_lock4.dat"),
            ],
            ignore_index=True,
        )
    elif days_irregular[day_irregular_index] == "20250321":
        shift_data_Sr = pd.concat(
            [
                open_shiftfile_Sr(path + days_irregular[day_irregular_index] + "/" + days_irregular[day_irregular_index] + "_clock_lock0.dat"),
                open_shiftfile_Sr(path + days_irregular[day_irregular_index] + "/" + days_irregular[day_irregular_index] + "_clock_lock1.dat"),
            ],
            ignore_index=True,
        )
else: ## For days (regular)
    shift_data_Sr = open_shiftfile_Sr(path + days[day_index] + "/" + days[day_index] + "_clock_lock0.dat")

maser_corrections = open_maser_correction(path + "daily maser offsets.csv")

if day_irregular_index != None:
    data_ErYb = open_ErYb_data(path + days_irregular[day_irregular_index] + "/" + days_irregular[day_irregular_index] + "_Deglitched_ErYb_only.dat") 
    shift_data_Al = open_shiftfile_Al(path + days_irregular[day_irregular_index] + "/" + days_irregular[day_irregular_index] + "_Alp_Freq_Shifts_ErYb.dat")
    shift_data_Yb = open_shiftfile_Yb(path + days_irregular[day_irregular_index] + "/YbI_1_rerun.txt")
else:
    data_ErYb = open_ErYb_data(path + days[day_index] + "/" + days[day_index] + "_Deglitched_ErYb_only.dat") 
    shift_data_Al = open_shiftfile_Al(path + days[day_index] + "/" + days[day_index] + "_Alp_Freq_Shifts_ErYb.dat")
    shift_data_Yb = open_shiftfile_Yb(path + days[day_index] + "/YbI_1_rerun.txt")


def compute_nuAl_ErYb(data):
    data["nuAl"] = -Decimal("105e6") + Decimal("560444") * (Decimal("1e9") + data["SDR:frep_ErYb"]) / Decimal(2) - data["fb_Al_ErYb"]
    data["nuAl"] = Decimal(4) * data["nuAl"]   

def compute_nuSr_ErYb(data):
    data["nuSi"] = -Decimal("105e6") + Decimal("388752") * (Decimal("1e9") + data["SDR:frep_ErYb"]) / Decimal(2) - Decimal("100e6")
    data["nuSr"] = (Decimal("1716882") / Decimal("777577")) * (data["nuSi"] - Decimal("216e6"))

def compute_nuYb_ErYb(data):
    data["nuYb"] = -Decimal("105e6") + Decimal("518237") * (Decimal("1e9") + data["SDR:frep_ErYb"]) / Decimal(2) - data["fb_Yb_ErYb"]
    data["nuYb"] = Decimal(2) * data["nuYb"] 

compute_nuAl_ErYb(data_ErYb)
compute_nuSr_ErYb(data_ErYb)
compute_nuYb_ErYb(data_ErYb) 

YbSrRatio2020 = Decimal("1.2075070393433378482") 
AlYbRatio2020 = Decimal("2.162887127516663703")
AlSrRatio2020 = Decimal("2.611701431781463025")
 
if day_irregular_index != None:
    correction_condition = days_irregular[day_irregular_index] == maser_corrections["date"]
else: 
    correction_condition = days[day_index] == maser_corrections["date"]
masercorrection = maser_corrections[correction_condition]["maser_offset"].apply(Decimal)

GR_shift_Al = Decimal("-8.114e-16") 
GR_shift_Yb = Decimal("-8.109e-16")
GR_shift_Sr = Decimal("10.660e-16")
GR_shift_sea_level = Decimal("-1798.501e-16")

total_correction_Yb = Decimal("1") + GR_shift_Yb + GR_shift_sea_level + masercorrection
total_correction_Sr = Decimal("1") + GR_shift_Sr + GR_shift_sea_level + masercorrection
total_correction_Al = Decimal("1") + GR_shift_Al + GR_shift_sea_level + masercorrection

#%%
### INITIAL SAMPLING RATES 
sampling_table_initial = [
    ["Al+", '{:0.2}'.format( np.mean(np.diff(shift_data_Al["MJD"]))*24*3600 ), np.percentile(np.diff(shift_data_Al["MJD"]), [0, 25, 50, 75, 100])*24*3600],
    ["Sr", '{:0.2}'.format( np.mean(np.diff(shift_data_Sr["MJD"]))*24*3600 ), np.percentile(np.diff(shift_data_Sr["MJD"]), [0, 25, 50, 75, 100])*24*3600],
    ["Yb", '{:0.2}'.format( np.mean(np.diff(shift_data_Yb["MJD"]))*24*3600 ), np.percentile(np.diff(shift_data_Yb["MJD"]), [0, 25, 50, 75, 100])*24*3600],
    ["ErYb", '{:0.2}'.format( np.mean(np.diff(data_ErYb["MJD"]))*24*3600 ), np.percentile(np.diff(data_ErYb["MJD"].astype(float)), [0, 25, 50, 75, 100])*24*3600],
]
if day_irregular_index != None:
    print("Initial/nominal sampling rates for ", days_irregular[day_irregular_index], "\n")
else: 
    print("Initial/nominal sampling rates for ", days[day_index], "\n")
print(tabulate(sampling_table_initial, headers=["Data", "Ave", "0   25   50   75    100"], tablefmt="grid"))


#%%
## SAMPLING RATES FOR FILTERED (COMB) DATA 
## filter comb data where one of the three clocks is missing 
##NOTE: this should be updated to filter according to clock pair (rather than triple) 
# comb_condition_AlSr = (~data_ErYb['nuAl'].isna() & ~data_ErYb['nuSr'].isna())
# comb_condition_YbSr = (~data_ErYb['nuSr'].isna() & ~data_ErYb['nuYb'].isna())
# comb_condition_AlYb = (~data_ErYb['nuAl'].isna() & ~data_ErYb['nuYb'].isna())
# comb_full_AlSr = data_ErYb[comb_condition_AlSr]
# comb_full_YbSr = data_ErYb[comb_condition_YbSr]
# comb_full_AlYb = data_ErYb[comb_condition_AlYb]

comb_condition = (~data_ErYb['nuAl'].isna() & ~data_ErYb['nuSr'].isna() & ~data_ErYb['nuYb'].isna())
comb_full = data_ErYb[comb_condition]

sampling_table_comb = [
    ["ErYb", '{:0.2}'.format( np.mean(np.diff(comb_full["MJD"]))*24*3600 ), np.percentile(np.diff(comb_full["MJD"].astype(float)), [0, 25, 50, 75, 100])*24*3600],
]
# sampling_table_comb = [
#     ["nuAl/nuSr", '{:0.2}'.format( np.mean(np.diff(comb_full_AlSr["MJD"]))*24*3600 ), np.percentile(np.diff(comb_full_AlSr["MJD"].astype(float)), [0, 25, 50, 75, 100])*24*3600],
#     ["nuYb/nuSr", '{:0.2}'.format( np.mean(np.diff(comb_full_YbSr["MJD"]))*24*3600 ), np.percentile(np.diff(comb_full_YbSr["MJD"].astype(float)), [0, 25, 50, 75, 100])*24*3600],
#     ["nuAl/nuYb", '{:0.2}'.format( np.mean(np.diff(comb_full_AlYb["MJD"]))*24*3600 ), np.percentile(np.diff(comb_full_AlYb["MJD"].astype(float)), [0, 25, 50, 75, 100])*24*3600]
# ]
if day_irregular_index != None:
    print("Sampling rates for non-missing data ", days_irregular[day_irregular_index], "\n")
else: 
    print("Sampling rates for non-missing data ", days[day_index], "\n")
print(tabulate(sampling_table_comb, headers=["Data", "Ave", "0   25   50   75    100"], tablefmt="grid"))


#%%
### SAMPLING RATES FOR FILTERED (SHIFT) DATA 
## filter NA data
al_cond = ~shift_data_Al['shift'].isna()
Al_non_na = shift_data_Al[al_cond]
Al = pd.Series(Al_non_na['MJD'])
## filter flagged data
good_condition_al = Al_non_na["IS_GOOD"] == 1
shift_data_Al_good = Al_non_na[good_condition_al].reset_index(drop=True, inplace = False)
# Al_good = pd.Series(shift_data_Al_good['MJD']) ##i don't think this are used/needed in kalman_investigation.ipynb

sr_cond = ~shift_data_Sr['shift'].isna()
Sr_non_na = shift_data_Sr[sr_cond]
Sr = pd.Series(Sr_non_na['MJD'])
good_condition_sr = Sr_non_na["IS_GOOD"] == 1
shift_data_Sr_good = Sr_non_na[good_condition_sr].reset_index(drop=True, inplace = False)
# Sr_good = pd.Series(shift_data_Sr_good['MJD'])

yb_cond = ~shift_data_Yb['shift'].isna()
Yb_non_na = shift_data_Yb[yb_cond]
Yb = pd.Series(Yb_non_na['MJD']) 
good_condition_yb = Yb_non_na["IS_GOOD"] == 1
shift_data_Yb_good = Yb_non_na[good_condition_yb].reset_index(drop=True, inplace = False)
# Yb_good = pd.Series(shift_data_Yb_good['MJD'])


sampling_table_shift_filtered = [
    ["Al+", '{:0.2}'.format( np.mean(np.diff(shift_data_Al_good["MJD"]))*24*3600 ), np.percentile(np.diff(shift_data_Al_good["MJD"]), [0, 25, 50, 75, 100])*24*3600],
    ["Sr", '{:0.2}'.format( np.mean(np.diff(shift_data_Sr_good["MJD"]))*24*3600 ), np.percentile(np.diff(shift_data_Sr_good["MJD"]), [0, 25, 50, 75, 100])*24*3600],
    ["Yb", '{:0.2}'.format( np.mean(np.diff(shift_data_Yb_good["MJD"]))*24*3600 ), np.percentile(np.diff(shift_data_Yb_good["MJD"]), [0, 25, 50, 75, 100])*24*3600],
]
if day_irregular_index != None:
    print("Sampling rates for filtered shift data on ", days_irregular[day_irregular_index], "\n")
else: 
    print("Sampling rates for filtered shift data on ", days[day_index], "\n")
print(tabulate(sampling_table_shift_filtered, headers=["Data", "Ave", "0   25   50   75    100"], tablefmt="grid"))





#%%
### FIND OVERLAPPING WINDOW
##TODO: generalize this to look for paired windows 

len_comb = len(comb_full['MJD']) 
len_Al = len(shift_data_Al_good['shift'])        
len_Sr = len(shift_data_Sr_good['shift'])        
len_Yb = len(shift_data_Yb_good['shift'])

def lb_extract(target, data):
    inx = 0
    stopper = 1
    while stopper == 1:
        if data.iloc[inx] < target:
            inx += 1
        else:
            return inx  

def ub_extract(target, data):
    inx = 1
    stopper = 1
    while stopper == 1:
        if data.iloc[len(data)-inx] > target:
            inx += 1
        else:
            return len(data)-inx 

# print("Comb start and end MJD: [", '{:0.11}'.format(comb_full['MJD'].iloc[0]), ', ', '{:0.11}'.format(comb_full['MJD'].iloc[len_comb-1]), ']')
# print("Al good shift start and end MJD: [", shift_data_Al_good['MJD'].iloc[0], ', ', shift_data_Al_good['MJD'].iloc[len_Al-1], ']')
# print("Sr good shift start and end MJD: [", shift_data_Sr_good['MJD'].iloc[0], ', ', shift_data_Sr_good['MJD'].iloc[len_Sr-1], ']')
# print("Yb good shift start and end MJD: [", shift_data_Yb_good['MJD'].iloc[0], ', ', shift_data_Yb_good['MJD'].iloc[len_Yb-1], ']')

starts = [comb_full['MJD'].iloc[0], shift_data_Al_good['MJD'].iloc[0], shift_data_Sr_good['MJD'].iloc[0], shift_data_Yb_good['MJD'].iloc[0]] 
ends = [comb_full['MJD'].iloc[len_comb-1], shift_data_Al_good['MJD'].iloc[len_Al-1], shift_data_Sr_good['MJD'].iloc[len_Sr-1], shift_data_Yb_good['MJD'].iloc[len_Yb-1]] 

last_start_time = max(starts)
first_end_time = min(ends)
print("Last start time: ", last_start_time)
print("First end time: ", first_end_time)

comb_start = ub_extract(target = last_start_time, data = comb_full['MJD'])  
comb_end = lb_extract(target = first_end_time, data = comb_full['MJD']) 

comb = pd.DataFrame()
comb["MJD"] = comb_full['MJD'].iloc[comb_start:comb_end] 
comb["nuAl"] = comb_full['nuAl'].iloc[comb_start:comb_end]
comb["nuSr"] = comb_full['nuSr'].iloc[comb_start:comb_end]
comb["nuYb"] = comb_full['nuYb'].iloc[comb_start:comb_end]
comb.reset_index(drop=True, inplace=True)

al_start = ub_extract(target = last_start_time, data = shift_data_Al_good["MJD"])
al_end = lb_extract(target = first_end_time, data = shift_data_Al_good["MJD"])  
shift_data_Al = shift_data_Al_good[al_start:al_end] 
shift_data_Al.reset_index(drop=True, inplace=True)

sr_start = ub_extract(target = last_start_time, data = shift_data_Sr_good["MJD"])
sr_end = lb_extract(target = first_end_time, data = shift_data_Sr_good["MJD"])  
shift_data_Sr = shift_data_Sr_good[sr_start:sr_end]
shift_data_Sr.reset_index(drop=True, inplace=True)

yb_start = ub_extract(target = last_start_time, data = shift_data_Yb_good["MJD"])
yb_end = lb_extract(target = first_end_time, data = shift_data_Yb_good["MJD"])  
shift_data_Yb = shift_data_Yb_good[yb_start:yb_end]
shift_data_Yb.reset_index(drop=True, inplace=True)

print("nuAl, nuSr, and nuYb start and end MJD: [", '{:0.11}'.format(comb["MJD"].iloc[0]), ', ', '{:0.11}'.format(comb["MJD"].iloc[len(comb["MJD"])-1]), ']')
print("Al good shift start and end MJD: [", shift_data_Al['MJD'].iloc[0], ', ', shift_data_Al['MJD'].iloc[len(shift_data_Al['MJD'])-1], ']')
print("Sr good shift start and end MJD: [", shift_data_Sr['MJD'].iloc[0], ', ', shift_data_Sr['MJD'].iloc[len(shift_data_Sr['MJD'])-1], ']')
print("Yb good shift start and end MJD: [", shift_data_Yb['MJD'].iloc[0], ', ', shift_data_Yb['MJD'].iloc[len(shift_data_Yb['MJD'])-1], ']')