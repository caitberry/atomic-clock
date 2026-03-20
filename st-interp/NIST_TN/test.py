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
day_irregular_index = 3 ##set to None if analyzing days, set to 0:5 if analyzing days_irregular

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
else: 
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
# print(total_correction_Yb)
# print(total_correction_Sr)
# print(total_correction_Al)

al_cond = ~shift_data_Al['shift'].isna()
Al_non_na = shift_data_Al[al_cond]
Al = pd.Series(Al_non_na['MJD'])

sr_cond = ~shift_data_Sr['shift'].isna()
Sr_non_na = shift_data_Sr[sr_cond]
Sr = pd.Series(Sr_non_na['MJD'])

yb_cond = ~shift_data_Yb['shift'].isna()
Yb_non_na = shift_data_Yb[yb_cond]
Yb = pd.Series(Yb_non_na['MJD']) 