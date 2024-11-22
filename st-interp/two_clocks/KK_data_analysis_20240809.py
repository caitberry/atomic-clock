#! /usr/local/bin/python
# %%
import pandas as pd
from decimal import Decimal, getcontext
import os
import matplotlib.pyplot as plt

# Set precision
getcontext().prec = 22

# Constants
masercorrection = Decimal("7.36631e-12")
fSr_BIPM = Decimal("429228004229872.99")
fYb_BIPM = Decimal("518295836590863.63")

total_correction_Yb = (
    Decimal("1") - Decimal("8.109e-16") - Decimal("1798.501e-16") - masercorrection
)
total_correction_Sr = (
    Decimal("1") + Decimal("10.660e-16") - Decimal("1798.501e-16") - masercorrection
)


print("INFO: total correction Yb: ", total_correction_Yb - Decimal(1))
print("INFO: total correction Sr: ", total_correction_Sr - Decimal(1))

ratio_BIPM = fYb_BIPM / fSr_BIPM
print(
    "INFO: BIPM - BACON: {:.3e}".format(ratio_BIPM - Decimal("1.2075070393433378482"))
)

# data diagonotics path
#data_diagonistics_path = "/Users/smt3/Documents/test-python"
##data_diagonistics_path = os.path.join("..", "data", "analysis_comparison", "KK")
#if not os.path.exists(data_diagonistics_path):
#    os.makedirs(data_diagonistics_path)

# %% Functions
def open_ErYb_data(data_path):
    # keys to read out as string
    key2read = ["fo_ErYb", "SDR:fb_Yb_ErYb", "SDR:frep_ErYb", "fb_Si_ErYb"]
    types = {key: str for key in key2read}
    types["MJD"] = float

    # # Read the CSV file
    data = pd.read_csv(
        data_path, header=1, delimiter="\t", dtype=types, engine="python"
    )
    # print(data.keys)
    # Convert the strings to Decimal for the given keys
    for k in key2read:
        data[k] = data[k].apply(Decimal)

    # reindex data
    data.index = range(len(data))

    return data[list(types.keys())]


def open_TiSa_data(data_path):
    # keys to read out as string
    key2read = ["fo_TiS", "fb_Yb_TiS", "SDR:frep_TiS", "fb_Si_TiS"]
    types = {key: str for key in key2read}
    types["MJD"] = float

    # Read the CSV file
    data = pd.read_csv(
        data_path, header=1, delimiter="\t", dtype=types, engine="python"
    )

    # Convert the strings to Decimal for the given keys
    for k in key2read:
        data[k] = data[k].apply(Decimal)

    # reindex data
    data.index = range(len(data))

    return data[list(types.keys())]


def open_shiftfile_Sr(datapath):
    data = pd.read_csv(
        datapath, header=21, delimiter="\t", dtype={1: str}, engine="python"
    )

    # Replace column names
    data.columns = ["MJD", "shift", "IS_GOOD"]

    # Change column type from float to bool
    data["IS_GOOD"] = data["IS_GOOD"].apply(lambda x: x == 1.0)

    # Remove rows with IS_GOOD column of 0
    data = data[data["IS_GOOD"] == True]

    # Convert the strings to Decimal for high precision
    data["shift"] = data["shift"].apply(Decimal)

    # reindex the data
    data.index = range(len(data))

    return data


def open_shiftfile_Yb(datapath):
    data = pd.read_csv(
        datapath, header=7, delimiter=r"\s{3,}", dtype={1: str}, engine="python"
    )

    # Replace column names
    data.columns = ["MJD", "shift", "IS_GOOD"]

    # Change column type from float to bool
    data["IS_GOOD"] = data["IS_GOOD"].apply(lambda x: x == 1.0)

    # Remove rows with IS_GOOD column of 0
    data = data[data["IS_GOOD"] == True]

    # Convert the strings to Decimal for high precision
    data["shift"] = data["shift"].apply(Decimal)

    # reindex the data
    data.index = range(len(data))

    return data


def compute_nuSr_TiSa(data):
    data["nuSi"] = (
        data["fo_TiS"]
        + Decimal("388944") * (Decimal("1e9") - data["SDR:frep_TiS"])
        - data["fb_Si_TiS"]
    )
    data["nuSr"] = (Decimal("1716882") / Decimal("777577")) * (
        data["nuSi"] / Decimal("2") - Decimal("248e6")
    )
    # remove rows with NaNs
    data.dropna(inplace=True)

    # reindex
    data.index = range(len(data))


def compute_nuSr_ErYb(data):
    data["nuSi"] = (
        -data["fo_ErYb"]
        + Decimal("388752") * (Decimal("1e9") + data["SDR:frep_ErYb"]) / Decimal(2)
        - data["fb_Si_ErYb"]
    )
    data["nuSr"] = (Decimal("1716882") / Decimal("777577")) * (
        data["nuSi"] - Decimal("216e6")
    )
    # remove rows with NaNs
    data.dropna(inplace=True)

    # reindex
    data.index = range(len(data))


def compute_nuYb_ErYb(data):
    data["nuYb"] = (
        -data["fo_ErYb"]
        + Decimal("518237") * (Decimal("1e9") + data["SDR:frep_ErYb"]) / Decimal(2)
        - data["SDR:fb_Yb_ErYb"]
    )
    # Remove rows with NaNs
    data.dropna(inplace=True)


def compute_nuYb_TiSa(data):
    data["nuYb"] = (
        data["fo_TiS"]
        + Decimal("259246") * (Decimal("1e9") - data["SDR:frep_TiS"])
        + data["fb_Yb_TiS"]
    )
    # Remove rows with NaNs
    data.dropna(inplace=True)


def get_ratio(dat_Yb, dat_Sr, data_comb):
    # cycle time for each measurement
    delta_mjd_Yb = float(dat_Yb.loc[1, "MJD"] - dat_Yb.loc[0, "MJD"])
    delta_mjd_Sr = float(dat_Sr.loc[1, "MJD"] - dat_Sr.loc[0, "MJD"])
    # delta_mjd_ErYb = data_comb.loc[1, "MJD"] - data_comb.loc[0, "MJD"]

    # compute the frequency ratio based on comb data's time stamp (mjd)
    nuYb = []
    nuSr = []
    shiftYb = []
    shiftSr = []
    common_mjd = []
    for ii in range(len(data_comb["MJD"])):
        mjd = data_comb.loc[ii, "MJD"]
        # first, find two closest Yb's mjd from the current comb's mjd. skip if there is no data before or after
        try:
            indx_Yb_after = dat_Yb[dat_Yb["MJD"] > mjd].index[0]
            indx_Yb_before = dat_Yb[dat_Yb["MJD"] < mjd].index[-1]

            # similar for Sr
            indx_Sr_after = dat_Sr[dat_Sr["MJD"] > mjd].index[0]
            indx_Sr_before = dat_Sr[dat_Sr["MJD"] < mjd].index[-1]
        except:
            continue

        # if two closest data points are too far apart, then skip
        if (
            dat_Yb.loc[indx_Yb_after, "MJD"] - dat_Yb.loc[indx_Yb_before, "MJD"]
        ) > delta_mjd_Yb * 1.2:
            continue
        if (
            dat_Sr.loc[indx_Sr_after, "MJD"] - dat_Sr.loc[indx_Sr_before, "MJD"]
        ) > delta_mjd_Sr * 1.2:
            continue

        # linear interpolate the shift data points onto the comb data
        a = Decimal(mjd) - Decimal(dat_Sr.loc[indx_Sr_before, "MJD"])
        b = Decimal(dat_Sr.loc[indx_Sr_after, "MJD"]) - Decimal(
            dat_Sr.loc[indx_Sr_before, "MJD"]
        )
        y1 = Decimal(dat_Sr.loc[indx_Sr_before, "shift"])
        y2 = Decimal(dat_Sr.loc[indx_Sr_after, "shift"])
        shiftSr.append(a / b * (y2 - y1) + y1)

        a = Decimal(mjd) - Decimal(dat_Yb.loc[indx_Yb_before, "MJD"])
        b = Decimal(dat_Yb.loc[indx_Yb_after, "MJD"]) - Decimal(
            dat_Yb.loc[indx_Yb_before, "MJD"]
        )
        y1 = Decimal(dat_Yb.loc[indx_Yb_before, "shift"])
        y2 = Decimal(dat_Yb.loc[indx_Yb_after, "shift"])
        shiftYb.append(a / b * (y2 - y1) + y1)

        # append to common mjd
        common_mjd.append(mjd)
        nuYb.append(data_comb.loc[ii, "nuYb"])
        nuSr.append(data_comb.loc[ii, "nuSr"])

    return (
        common_mjd,
        nuYb,
        nuSr,
        shiftYb,
        shiftSr,
    )


# Nick' decimal average
def meanDecimal(arr):
    if (type(arr) == list) or (type(arr) == pd.core.series.Series):
        arr = [x for x in arr if not x.is_nan()]  # get rid of NaNs
        s = Decimal("0")
        for i in arr:
            s = s + i

        length = len(arr) if len(arr) > 0 else 1
        return s / Decimal(length)
    else:
        return arr


# %% Main

# data loading
data_ErYb = pd.concat(
    [
        open_ErYb_data(
            os.path.join("..", "smt3", "Documents", "test-python", "20240716_Deglitched_ErYb_only.dat")
        ),
    ],
    ignore_index=True,
)

data_TiSa = pd.concat(
    [
        open_TiSa_data(
            os.path.join("..", "smt3", "Documents", "test-python", "20240716_Deglitched_TiS_only.dat")
        )
    ],
    ignore_index=True,
)

shift_data_Sr = pd.concat(
    [
        open_shiftfile_Sr(
            os.path.join("..", "smt3", "Documents", "test-python",  "20240716_Sr_Freq_Shifts_1.dat")
        ),
        open_shiftfile_Sr(
            os.path.join("..", "smt3", "Documents", "test-python", "20240716_Sr_Freq_Shifts_2.dat")
        ),
        open_shiftfile_Sr(
            os.path.join("..", "smt3", "Documents", "test-python",  "20240716_Sr_Freq_Shifts_3.dat")
        ),
        open_shiftfile_Sr(
            os.path.join("..", "smt3", "Documents", "test-python",  "20240716_Sr_Freq_Shifts_4.dat")
        ),
    ],
    ignore_index=True,
)
shift_data_Yb = pd.concat(
    [
        open_shiftfile_Yb(
            os.path.join("..", "smt3", "Documents", "test-python",   "20240716_Yb_Freq_Shifts.txt")
        )
    ],
    ignore_index=True,
)

# data processing
compute_nuSr_ErYb(data_ErYb)
compute_nuSr_TiSa(data_TiSa)
compute_nuYb_ErYb(data_ErYb)
compute_nuYb_TiSa(data_TiSa)

# Calculate optical frequenceis
common_mjd_ErYb, nuYb_ErYb, nuSr_ErYb, shiftYb_ErYb, shiftSr_ErYb = get_ratio(
    shift_data_Yb, shift_data_Sr, data_ErYb
)
common_mjd_TiSa, nuYb_TiSa, nuSr_TiSa, shiftYb_TiSa, shiftSr_TiSa = get_ratio(
    shift_data_Yb, shift_data_Sr, data_TiSa
)


# %% Outputs
frequency_Yb_ErYb = [
    (Decimal(2) * nuYb_ErYb[i] + shiftYb_ErYb[i]) * total_correction_Yb
    for i in range(len(nuYb_ErYb))
]
frequency_Sr_ErYb = [
    (nuSr_ErYb[i] + shiftSr_ErYb[i]) * total_correction_Sr
    for i in range(len(nuSr_ErYb))
]
frequency_ratio_ErYb = [
    frequency_Yb_ErYb[i] / frequency_Sr_ErYb[i] for i in range(len(frequency_Yb_ErYb))
]

frequency_Yb_TiSa = [
    (Decimal(2) * nuYb_TiSa[i] + shiftYb_TiSa[i]) * total_correction_Yb
    for i in range(len(nuYb_TiSa))
]
frequency_Sr_TiSa = [
    (nuSr_TiSa[i] + shiftSr_TiSa[i]) * total_correction_Sr
    for i in range(len(nuSr_TiSa))
]
frequency_ratio_TiSa = [
    frequency_Yb_TiSa[i] / frequency_Sr_TiSa[i] for i in range(len(frequency_Yb_TiSa))
]

# Calculate ratio differences
ratio_difference_ErYb = [x - ratio_BIPM for x in frequency_ratio_ErYb]
mean_diff_ErYb = meanDecimal(ratio_difference_ErYb)

ratio_difference_TiSa = [x - ratio_BIPM for x in frequency_ratio_TiSa]
mean_diff_TiSa = meanDecimal(ratio_difference_TiSa)

# print mean differences with 3 significant digits
print("INFO: ratio difference ErYb: ", mean_diff_ErYb)
print("INFO: ratio difference TiSa: ", mean_diff_TiSa)

# Calculate mean frequencies each
mean_freq_Sr_ErYb = meanDecimal(frequency_Sr_ErYb)
mean_freq_Yb_ErYb = meanDecimal(frequency_Yb_ErYb)
mean_freq_Sr_TiSa = meanDecimal(frequency_Sr_TiSa)
mean_freq_Yb_TiSa = meanDecimal(frequency_Yb_TiSa)

# print deviation from BIPM
print(
    "INFO: Yb frequency deviation from BIPM, ErYb: ",
    (mean_freq_Yb_ErYb - fYb_BIPM) / fYb_BIPM,
)
print(
    "INFO: Sr frequency deviation from BIPM, ErYb: ",
    (mean_freq_Sr_ErYb - fSr_BIPM) / fSr_BIPM,
)
print(
    "INFO: Yb frequency deviation from BIPM, TiSa: ",
    (mean_freq_Yb_TiSa - fYb_BIPM) / fYb_BIPM,
)
print(
    "INFO: Sr frequency deviation from BIPM, TiSa: ",
    (mean_freq_Sr_TiSa - fSr_BIPM) / fSr_BIPM,
)


fig, axes = plt.subplots(2, 1, figsize=(10, 10), dpi=200)
axes[0].plot(common_mjd_ErYb, ratio_difference_ErYb, ".", label="ErYb")
axes[0].set_xlabel("MJD")
axes[0].set_ylabel("ratio offset from BIPM")
axes[0].legend()
axes[0].grid()
axes[1].plot(common_mjd_TiSa, ratio_difference_TiSa, ".", label="TiSa")
axes[1].set_xlabel("MJD")
axes[1].set_ylabel("ratio offset from BIPM")
axes[1].legend()
axes[1].grid()
fig.tight_layout()
# plt.show()


# %% save some data for diagonotics
data_interpolated_ErYb = pd.DataFrame({})
data_interpolated_ErYb["MJD"] = common_mjd_ErYb
data_interpolated_ErYb["nuYb"] = [nuYb_ErYb[ii] * total_correction_Yb for ii in range(len(nuYb_ErYb))]
data_interpolated_ErYb["nuSr"] = [nuSr_ErYb[ii] * total_correction_Yb for ii in range(len(nuSr_ErYb))]
data_interpolated_ErYb["shiftYb"] = shiftYb_ErYb
data_interpolated_ErYb["shiftSr"] = shiftSr_ErYb
data_interpolated_ErYb["ratio"] = frequency_ratio_ErYb

data_interpolated_TiSa = pd.DataFrame({})
data_interpolated_TiSa["MJD"] = common_mjd_TiSa
data_interpolated_TiSa["nuYb"] = [nuYb_TiSa[ii] * total_correction_Yb for ii in range(len(nuYb_TiSa))]
data_interpolated_TiSa["nuSr"] = [nuSr_TiSa[ii] * total_correction_Yb for ii in range(len(nuSr_TiSa))]
data_interpolated_TiSa["shiftYb"] = shiftYb_TiSa
data_interpolated_TiSa["shiftSr"] = shiftSr_TiSa
data_interpolated_TiSa["ratio"] = frequency_ratio_TiSa

data_raw = pd.DataFrame({})
data_raw["MJD_ErYb"] = data_ErYb["MJD"]
data_raw["frep"] = data_ErYb["SDR:frep_ErYb"] * total_correction_Yb
data_raw["nuYb"] = data_ErYb["nuYb"] * total_correction_Yb
data_raw["nuSr"] = data_ErYb["nuSr"] * total_correction_Yb

data_raw["MJD_TiSa"] = data_TiSa["MJD"]
data_raw["frep"] = data_TiSa["SDR:frep_TiS"] * total_correction_Yb
data_raw["nuYb"] = data_TiSa["nuYb"] * total_correction_Yb
data_raw["nuSr"] = data_TiSa["nuSr"] * total_correction_Yb

data_raw["MJD_Yb"] = shift_data_Yb["MJD"]
data_raw["shiftYb"] = shift_data_Yb["shift"]

data_raw["MJD_Sr"] = shift_data_Sr["MJD"]
data_raw["shiftSr"] = shift_data_Sr["shift"]

# save the dataframes
data_interpolated_ErYb.to_csv(
    os.path.join(data_diagonistics_path, "data_interpolated_ErYb.csv"), index=False
)
data_interpolated_TiSa.to_csv(
    os.path.join(data_diagonistics_path, "data_interpolated_TiSa.csv"), index=False
)
data_raw.to_csv(os.path.join(data_diagonistics_path, "data_raw.csv"), index=False)
print("INFO: data saved to ", data_diagonistics_path)

# %%