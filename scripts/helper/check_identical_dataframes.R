


dataframe2check_01 <- data01
dataframe2check_02 <- data01_old

# Assuming data01 and data01_old are your data frames

# Initialize a vector to store names of columns that differ
different_columns <- c()

# Loop through the columns of data01
for (col_name in names(dataframe2check_01)) {
  # Check if the column exists in data01_old
  if (col_name %in% names(dataframe2check_02)) {
    # Compare the columns using all.equal for numerical values with a tolerance
    if (!isTRUE(all.equal(dataframe2check_01[[col_name]], dataframe2check_02[[col_name]], tolerance = 1e-8))) {
      # If not all.equal, add to the list of different columns
      different_columns <- c(different_columns, col_name)
    }
  } else {
    # If the column doesn't exist in data01_old, note it as different
    different_columns <- c(different_columns, col_name)
  }
}

# Print the names of columns that have differences
if (length(different_columns) > 0) {
  cat("Columns with differences:", paste(different_columns, collapse = ", "), "\n")
} else {
  cat("No differences found in columns.\n")
}

dataframe2check_01$auc[1:5]
dataframe2check_02$auc[1:5]
