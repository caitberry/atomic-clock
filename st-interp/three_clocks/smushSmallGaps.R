smushSmallGaps = function(dat) {
  # Find indices of missing (NA) values in the data.
  missing.indices <- which(is.na(dat))
  # Calculate the differences between consecutive missing indices.
  differenced <- diff(missing.indices)
  # Identify where the gap between consecutive missing indices is greater than 1.
  jumps <- which(differenced > 1)
 
  # Make a copy of the original data.
  new.dat <- dat
  # Replace NA values with Inf (infinity) to mark them clearly.
  new.dat[is.na(dat)] <- Inf
  # Find consecutive runs of values and their lengths in the modified data.
  runs <- rle(new.dat) #run length encoding 
  # Extract the lengths of runs of Inf values (representing the missing periods).
  NA_period_length <- runs$lengths[runs$values == Inf]
  # Initialize a vector to store indices of missing values that will be removed.
  remove.indices <- c()
  # Initialize a counter to track the starting index of each missing period.
  start <- 0
  # Loop through each period of missing values (represented by Inf).
  for (i in 1:length(NA_period_length)) {
    # If the length of the missing period is less than 100, mark it for removal.
    if (NA_period_length[i] < 100) {
      remove.indices <- append(remove.indices, missing.indices[(1 + start):sum(NA_period_length[1:i])])  ##replace w/ imputing 
      # Update the start index for the next period.
      start <- sum(NA_period_length[1:i])
    }
    # If the length of the missing period is greater than 100, just update the start index.
    if (NA_period_length[i] > 100) {
      start <- sum(NA_period_length[1:i])
    }
  }
  # Remove the small gaps (NA values with length < 100).
  new.dat <- dat[-c(remove.indices)]
  # Return the cleaned-up data.
  return(new.dat)
}
 
### Summary of what the function does:
#This function processes a vector (`dat`) by identifying periods of missing values 
#(represented as `NA`). It removes "small gaps" of missing values that are shorter 
#than 100 entries, while retaining longer gaps. It uses a combination of finding runs 
#of `NA` values, measuring their lengths, and excluding short runs. The cleaned data 
#is returned without the smaller gaps.