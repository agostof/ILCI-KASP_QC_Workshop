# 2024 ILCI Project Meeting
# Workshop: Better breeding operations with QC KASP markers
# Activity #3: Identifying F1 progeny using molecular markers
# Lead: Clara Cruet-Burgos

# SET your working directory to where you downloaded the data!
#setwd("~/Desktop/KASP_QC_Workshop 2/")

# Looking at your directory
getwd()

# Cleaning your environment
rm(list = ls())

## Loading functions that we will use
#Function to filter_markers
filter_markers <- function(results, Parent1, Parent2, CrossID) {
  # Subsetting all data for your selected cross
  verify_cross <- subset(results, grepl(CrossID, results$Crosses))
  
  # Selecting segregating markers
  verify_cross_seg <- verify_cross[1]  # Initialize the dataframe to store the result
  for (col_name in names(verify_cross)[4:ncol(verify_cross)]) {  # Looping over columns
    alleles <- length(unique(verify_cross[[col_name]]))  # Verify if they are segregating
    
    if (alleles > 1) {  # If they are segregating, keep the column
      # Keep the column if there are more than one unique alleles
      verify_cross_seg[[col_name]] <- verify_cross[[col_name]]  # Add it to new data
    }
  }
  
  # Subsetting parents and binding on top
  parent1_data <- results[results$Sample == Parent1, ]  # Extracting parent information
  parent2_data <- results[results$Sample == Parent2, ]  # Extracting parent information
  parental_data <- rbind(parent1_data, parent2_data)  # Merging parent information
  
  # Filter markers that were segregating in progeny
  parental_data_seg <- parental_data[names(parental_data) %in% colnames(verify_cross_seg)]
  
  # Data for decision, binding parent and crosses
  segregating_markers <- rbind(parental_data_seg, verify_cross_seg)
  
  return(segregating_markers)
}

#### Function to remove reference alleles
filter_reference_alleles <- function(segregating_markers, Parent1, Parent2, columns_to_remove) {
  # Remove specified columns from segregating_markers
  segregating_markers <- segregating_markers[, !names(segregating_markers) %in% columns_to_remove, drop = FALSE]
  
  # Set row names
  row.names(segregating_markers) <- segregating_markers$Sample
  segregating_markers <- segregating_markers[-1]
  
  # Define the reference rows (row names)
  reference_rows <- c(Parent1, Parent2)
  
  # Get the reference data
  reference_data <- cbind(segregating_markers[reference_rows, ])
  rownames(reference_data) <- reference_rows
  colnames(reference_data) <- names(segregating_markers) 
  
  # Create a logical condition to filter the dataframe
  condition <- apply(segregating_markers, 1, function(row) !any(apply(reference_data, 1, function(ref_row) all(row == ref_row))))
  
  # Filter the dataframe based on the condition
  filtered_results_new <- segregating_markers[condition, , drop = FALSE]
  
  # Combine reference_data and filtered_results_new
  candidate_f1 <- rbind(reference_data, filtered_results_new)
  
  return(candidate_f1)
}

# Loading the result data and doing some edits
results <- read.csv("ResultsGrid.csv", skip = 6)  # Skip the first six lines
names(results) <- results[1, ]  # Use the first row as header or column names
names(results)[1] <- 'Sample'  # Edit the first column name
results <- results[-1, ]  # Remove the first row
results$Crosses <- gsub("_.*", "", results$Sample)  # Create a new column with cross names
View(results)

# Verify your cross
Parent1 <- "PVT-305" # parental line 1
Parent2 <- "BTP-810" # parental line 2
CrossID <- "PVT-BTP-F1" # parental line 3

# Running filtering markers function to identify segregating markers
segregating_markers <- filter_markers(results, Parent1, Parent2, CrossID)
View(segregating_markers)

# Remove additional markers that might not be informative
columns_to_remove <- vector() # creating a filtering vector
columns_to_remove <- c("snpPV00087","snpPV00088", "snpPV00099") # markers to remove

# Running filtering reference alleles function to identify recombinants
candidate_f1 <- filter_reference_alleles(segregating_markers, Parent1, Parent2, columns_to_remove)
View(candidate_f1)
