# Data format for Google Drive 2024 clock comparison data folder
## owner: nicknardelli@gmail.com 

Note, I have copied the contents of the Google Drive folder and deleted all 2024 subfolders and renamed this master folder "2025 clock comparison data" to 
reflect that our interest is in the most recent 2025 data. From this copy of the Drive folder I have also removed the "Manuals_Menlo_portable" and "BRAN loopback 2025Jan". The "Red shifts.docx" and "daily maser offsets.xlsx" files contain information for the corrections component of the code. I changed the format of "daily maser offsets.xlsx" so that the date column is formatted as yyyy-mm-dd, renamed the column maser offset to "maser_offset" and formatted the data in this column so that it is understood as a scientific number with five decimal places. 

In addition to the data, there is a folder called "ratios" from a Wednesday, April 2, 2025 email from Nick to Amanda and Cait. The contents of this folder are used
to identify which dates are of interest in the final analysis. The dates and ratios are summarized in the following table. 

| Date | Clock 1 | Clock 2 |
| --- | --- | --- | 
| 20250227 | Al+ | Sr |
| 20250228 | Al+ | Sr |
| 20250304 | Al+ | Sr |
| 20250306 | Al+ | Sr |
| 20250307 | Al+ | Sr |
| 20250313 | Al+ | Sr |
| 20250318 | Al+ | Sr |
| 20250320 | Al+ | Sr |
| 20250321 | Al+ | Sr |
| 20250227 | Al+ | Yb |
| 20250228 | Al+ | Yb |
| 20250304 | Al+ | Yb |
| 20250306 | Al+ | Yb |
| 20250307 | Al+ | Yb |
| 20250313 | Al+ | Yb |
| 20250318 | Al+ | Yb |
| 20250320 | Al+ | Yb |
| 20250321 | Al+ | Yb |
| 20250116 | Yb | Sr |
| 20250124 | Yb | Sr |
| 20250204 | Yb | Sr |
| 20250206 | Yb | Sr |
| 20250227 | Yb | Sr |
| 20250228 | Yb | Sr |
| 20250304 | Yb | Sr |
| 20250306 | Yb | Sr |
| 20250307 | Yb | Sr |
| 20250313 | Yb | Sr |
| 20250318 | Yb | Sr |
| 20250320 | Yb | Sr |
| 20250321 | Yb | Sr |


## Original format of data 

All subfolders except the "YbReRunShifts" follow the date convention YYYYMMDD. 

Within each date folder there are:
    - .dat files corresponding to Al+, Sr clock lock, and deglitched comb 
    - .txt file of Yb data 
    - subfolder of various names, unsure of the purpose of these files and how they relate to the .dat and .txt files mentioned above 
        * shifts
        * shift_tables 
        * shift tables
        * systematics 
        * error budgets 
        * error tables
        * error budget 
The folder 20250313 also contains a subfolder called "Before filtering autolocker issue" containing 
    - a folder called "20250313_after_filtering"
    - a folder called "20250313_before_filtering"
    - two comb .dat files 
The folder YbReRunShifts contains another folder called "YbReEvaluation_250403" which contains dated subfolders. The names of these subfolders follow the date convention
YYMMDD. These dated subfolders match all of the data subfolders within the main folder except there is no match for the date 2025-03-14. 
These dated folders each contain two files:
    - YbI_1.txt or YbI_2.txt
    - YbI_StaticShift1.txt or YbI_StaticShift2.txt 
For the Sr data, all dates besides the following have only clock_lock0 data files:
    - 20250206 has lock0 and lock2; 
    - 20250228 has lock1 only; 
    - 20250306 has lock2 only; 
    - 20250313 has lock0, lock1, lock2, and lock3; 
    - 20250320 has lock0, lock2, and lock4; 
    - 20250321 has lock0 and lock1. 

## Edited format of data 
All files corresponding to the same date have been put in the same folder following the naming convention YYYYMMDD. There is no more YbReRunShifts folder. The 
re-run .txt files for Yb have been renamed "YbI_1_rerun.txt" or "YbI_2_rerun.txt". I have deleted data from 20240314 since this date does not appear in Nick's ratio data files (see table above).

All subfolders that were originally named some variation of shift or error tables or systematics or error budget have all been renamed "systematics". The folders contain only .csv or .json files. I am stil uncertain how these files are related the the other data files in the same date parent folder. 

Renamed all clock_lock#X.dat files to remove # sign. Renamed all Al+_Freq_Shifts_ErYb.dat files to replace + with p. 

The two YbI_2_rerun.txt files occured for dates 2025-02-04 and 2025-03-18 and have been renamed as YbI_1_rerun.txt. 

The hidden comb data file in the "Before filtering autolocker issue" subfolder of 2025-03-13 is formatted differently than the other comb data files (last two columns switched). 

For each folder of data (besides the 2025-03-13 "Before filtering autolocker issue" folder) I have verified that 
    - Each Deglitched_ErYb_only data file has two header lines to ignore and contains the columns: MJD, timer, SDR:frep_ErYb, fo_ErYb, fb_Si_ErYb, fb_Yb_ErYb, fb_Al_ErYb;
    - Each Alp_Freq_Shifts_ErYb data file has only 30 header lines to ignore and contains the columns: t_mjd, f_shift [Hz], IS_GOOD; 
    - Each clock_lock data file has only 24 header lines to ignore and contains the columns: t_mjd, f_shift [Hz], IS_GOOD; 
    - Each YbI_1_rerun data file has only 8 header lines to ignore and contains the columns:  t_mjd, f_shift [Hz], IS_GOOD. 
For the following dates this required edits to the original files to match the format above: 2025-02-27, 2025-02-28, 2025-03-04, 2025-03-06, 2025-03-07, 2025-03-13, 2025-03-18, 2025-03-20, 2020-03-21.  

### Initial analysis 

Restrict our calculations to dealing with only the following files:
    - Alp_Freq_Shifts_ErYb.dat 
    - clock_lock0.dat  
    - Deglitched_ErYb_only.dat 
    - YbI_1_rerun.txt 

Restrict initial analysis to all dates besides: 20250206, 20250228, 20250306, 20250313, 20250320, 20250321 because they require special data loading fundtions using pd.concat. For example: 

 shift_data_Sr = pd.concat(
    [
        open_shiftfile_Sr(path + "20240716_Sr_Freq_Shifts_1.dat"),
        open_shiftfile_Sr(path + "20240716_Sr_Freq_Shifts_2.dat"),
        open_shiftfile_Sr(path + "20240716_Sr_Freq_Shifts_3.dat"),
        open_shiftfile_Sr(path + "20240716_Sr_Freq_Shifts_4.dat"),
    ],
    ignore_index=True,
)

