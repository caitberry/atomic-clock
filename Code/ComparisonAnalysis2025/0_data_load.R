###### Load in Ratio Data #######

folderLocation="/home/aak3/NIST/atomic-clock/"

# List the folders containing .txt files
folders <- c(paste(folderLocation,"Data/ClockComparison2025/AlSr",sep=""), 
             paste(folderLocation,"Data/ClockComparison2025/AlYb",sep=""), 
             paste(folderLocation,"Data/ClockComparison2025/YbSr",sep=""))

# Initialize a list to store data from each file
data_list <- list()

# Loop through each folder and load .txt files
for (folder in folders) {
  # Get all .txt files in the folder
  file_paths <- list.files(path = folder, pattern = "\\.txt$", full.names = TRUE)
  
  # Read each file and store it in the list
  for (file in file_paths) {
    data <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)  # Adjust sep if needed
    
    # Extract the file name
    file_name <- basename(file)
    
    # Extract the date (assuming format YYYYMMDD at start)
    date_string <- regmatches(file_name, regexpr("^\\d{8}", file_name))
    data$date <- as.Date(date_string, format = "%Y%m%d")
    
    # Extract the ratio (text between the dash and ".txt")
    ratio_string <- sub("^\\d{8}-([^.]+)\\.txt$", "\\1", file_name)
    data$ratio <- ratio_string
    
    # Store data with file path as key
    data_list[[file]] <- data
  }
}
# 
# # Now, data_list contains all the loaded datasets
# # You can access individual datasets using their file paths as names
# 
# # Example: View the first few rows of a specific dataset
# # head(data_list[["folder1/file1.txt"]])
# 
# # If you prefer separate named objects instead of a list, you can assign them dynamically:
# for (name in names(data_list)) {
#   cleaned_name <- gsub("[^a-zA-Z0-9]", "_", basename(name))  # Replace non-alphanumeric chars with "_"
#   
#   # Ensure the name starts with a letter by prefixing with "data_" if necessary
#   if (grepl("^[0-9]", cleaned_name)) {
#     cleaned_name <- paste0("data_", cleaned_name)
#   }
#   
#   assign(cleaned_name, data_list[[name]], envir = .GlobalEnv)
# }
# 
# rm(data)
