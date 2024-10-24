# Install 'maftools' if not already installed
if (!requireNamespace("maftools", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("maftools")
}

setwd("/home/Ankita.Singh/Desktop/NEOANTIGEN/neoantigent/")

# Load the maftools library
library(maftools)
library(dplyr)

# Read the MAF file into 'all_df'
all_df <- read.maf(maf = "/home/Ankita.Singh/Desktop/NEOANTIGEN/Data/Others/mc3.v0.2.8.PUBLIC.maf.gz")

# Extract the gene summary from the MAF object
gene_summary <- getGeneSummary(all_df)

# Print the column names to ensure we're using the correct names
print(colnames(gene_summary))

# Define a function to calculate the difference and log for each gene
calculate_difference_and_log <- function(total, mutated) {
  # Subtract 'MutatedSamples' from 'total'
  if (is.na(total) | is.na(mutated)) {
    return(NA)
  }
  
  diff <- (total / mutated)* log(mutated+1)
  
  return(diff)
}
  # # Compute the natural log of the difference, return NA if diff <= 0
  # if (diff > 0) {
  #   log_diff <- log(diff)
  # } else {
  #   log_diff <- NA
  # }
  
#   return(log_diff)
# }
  

# Once we identify the correct column names, use them in this step
# Assuming the columns are correctly named (adjust after inspecting the output of colnames)

gene_summary <- gene_summary %>%
  rowwise() %>%
  mutate(log_difference = calculate_difference_and_log(total, MutatedSamples))  # Use correct column names here

# Display the genes and the computed log differences
gene_summary_result <- gene_summary %>% select(Hugo_Symbol, total, MutatedSamples, log_difference) %>% data.frame()
gene_summary_result <- gene_summary_result[order(gene_summary_result$log_difference, decreasing = T),]

# View the first few rows of the result
print(head(gene_summary_result))
# Save the entire gene summary to a CSV file
write.csv(gene_summary, "../Results/gene_summary_full.csv", row.names = FALSE)
