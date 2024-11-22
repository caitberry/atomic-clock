import pandas as pd
from decimal import Decimal
import matplotlib.pyplot as plt
import numpy as np
import csv
 
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
    data["IS_GOOD"] = data["IS_GOOD"].apply(lambda x: x == 1.0)
 
    # Put NaN in data["shift"] where data["IS_GOOD"] is 0
    data.loc[~data["IS_GOOD"], "shift"] = np.nan
 
    # Change column type to float
    data["shift"] = data["shift"].apply(float)
    data["MJD"] = data["MJD"].apply(float)
 
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
    data["MJD"] = data["MJD"].apply(float)
 
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
 
# load Yb shift data - note .txt not .dat may be causing issue w/ column type
shift_data_Yb = open_shiftfile_Yb(path + "20240716_Yb_Freq_Shifts.txt")
 
shift_data_Yb["shift"] = [Decimal(i) for i in shift_data_Yb["shift"]]
 
 
 
 
 
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
# and change type to high-precision Decimal
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
shift_data_Sr["shift"] = [Decimal(i) for i in shifts1s]
 


shiftSr = shift_data_Sr
shiftYb = shift_data_Yb

#print(shiftSr)

shiftSr.to_csv('Sr_shift_data.csv', index=False, na_rep= "NaN")
shiftYb.to_csv('Yb_shift_data.csv', index=False, na_rep= "NaN")
#data_ErYb.to_csv('ErYb_data.csv', index=False)

# np.savetxt("Sr_shift_data.csv", shiftSr, delimiter=",")
 
"""with open('Sr_shift_data.csv', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(shiftSr)

with open('Yb_shift_data.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(shiftYb)

with open('ErYb_data.csv', newline='\n') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(data_ErYb)"""


