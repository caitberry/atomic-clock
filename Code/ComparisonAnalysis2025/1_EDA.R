######## Explore Data ########
library(tidyverse)
library(magrittr)


## Put all ratio data for each clock into a single data frame with columns:
##  V1 == Modified Julian Date
##  V2 == ratio measurement
##  source_file == # day in order of ratio folder

library(dplyr)  # Load dplyr for binding rows easily

# Create a list to store the combined datasets by folder
combined_data <- list()

for (folder in folders) {
  # Filter out files that belong to the current folder
  file_paths <- names(data_list)[grepl(paste0("^", folder, "/"), names(data_list))]
  
  # Extract numbers from filenames for sorting
  file_numbers <- as.numeric(gsub("\\D", "", basename(file_paths)))  # Remove non-numeric characters
  
  # Sort files by extracted number
  sorted_files <- file_paths[order(file_numbers)]
  
  # Bind datasets row-wise
  combined_data[[folder]] <- bind_rows(lapply(sorted_files, function(f) data_list[[f]]), .id = "source_file")
}

# The combined_data list now contains one dataset per folder
# Example: View the first few rows of a merged dataset
head(combined_data[[paste(folderLocation,"Data/ClockComparison2025/AlSr",sep="")]])

AlSr_df <- combined_data[[paste(folderLocation,"Data/ClockComparison2025/AlSr",sep="")]]
YbSr_df <- combined_data[[paste(folderLocation,"Data/ClockComparison2025/YbSr",sep="")]]
AlYb_df <- combined_data[[paste(folderLocation,"Data/ClockComparison2025/AlYb",sep="")]]

# AlSr_df$source_file %<>% as.numeric
# YbSr_df$source_file %<>% as.numeric
# AlYb_df$source_file %<>% as.numeric

# par(mfrow = c(3,1))
# plot(AlYb_df$V1, AlYb_df$V2, pch = 19, xlab = "MJD", ylab = "AlYb Ratios")
# plot(AlSr_df$V1, AlSr_df$V2, pch = 19, xlab = "MJD", ylab = "AlSr Ratios")
# #changed the last one so that would only plot the same days as above clock ratios
# plot(YbSr_df$V1[as.numeric(YbSr_df$source_file) > 4], YbSr_df$V2[as.numeric(YbSr_df$source_file) > 4], pch = 19, xlab = "MJD", ylab = "YbSr Ratios")
# 
# 
