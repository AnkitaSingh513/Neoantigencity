# Install necessary packages if not already installed
if (!requireNamespace("maftools", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("maftools")
}

if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}

# Load required libraries
library(maftools)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)

setwd("/home/Ankita.Singh/Desktop/NEOTANTIGEN/")

# Step 1: Load the MAF file (mutation data)
maf_file <- "/home/Ankita.Singh/Desktop/NEOTANTIGEN/Data/Others/mc3.v0.2.8.PUBLIC.maf.gz"
maf_data <- read.maf(maf = maf_file)
##################################################

# Step 2: Load the GDC sample file
gdc_file <- "/home/Ankita.Singh/Desktop/NEOTANTIGEN/neoantigent/GDC-PANCAN.basic_phenotype.tsv"
gdc_data <- read_tsv(gdc_file)  # Use read_tsv for tab-delimited files
##################################################

# Step 3: Truncate 'sample' and 'Tumor_Sample_Barcode' to the first 12 characters
gdc_data <- gdc_data %>%
  mutate(sample_trimmed = substr(sample, 1, 12))  # Extract the first 12 characters from the GDC data
#Ragh is modifying
gdc_data$sample_trimmed <- gdc_data$sample
maf_data@data$Tumor_Sample_Barcode <- as.character(as.vector(maf_data@data$Tumor_Sample_Barcode))
temp <- strsplit(maf_data@data$Tumor_Sample_Barcode,"-")
tumor_samples <- paste(unlist(map(temp,c(1))),unlist(map(temp,c(2))),unlist(map(temp,c(3))),unlist(map(temp,c(4))),sep = "-")
maf_data@data$Tumor_Short_Barcode <- tumor_samples

mutation_data <- maf_data@data %>%
  select(Tumor_Short_Barcode, Hugo_Symbol, Variant_Classification, dbSNP_RS, Chromosome) %>%
  #mutate(Tumor_Sample_Barcode_trimmed = substr(Tumor_Sample_Barcode, 1, 12)) %>%
  filter(!is.na(Tumor_Short_Barcode), !is.na(Hugo_Symbol), !is.na(Variant_Classification))

# Step 4: Merge the mutation data with GDC sample data
merged_data <- mutation_data %>%
  left_join(gdc_data, by = c("Tumor_Short_Barcode" = "sample_trimmed"))  %>% #, relationship = "many-to-many")
  filter(!is.na(project_id))

# Step 5: Calculate global mutation frequencies
global_freq <- merged_data %>%
  group_by(Variant_Classification) %>%
  summarise(global_count = n()) #%>%
  #mutate(global_weight = global_count / sum(global_count))  # Normalize to get global weight

## Step 6: Group by necessary columns and calculate mutation counts
mutation_counts <- merged_data %>%
  group_by(Tumor_Short_Barcode, Hugo_Symbol, Variant_Classification, dbSNP_RS, Chromosome, project_id) %>%
  summarise(mutation_count = n(), .groups = 'drop')

## Step 7: Merge mutation counts with global frequencies to calculate weighted score
#mutation_scores <- mutation_counts %>%
#  left_join(global_freq, by = "Variant_Classification") %>%
#  mutate(weighted_score = mutation_count * global_weight)

# Step 8: Create a wide-format matrix including all specified mutation types as columns
# Include all the required columns and mutation types
final_matrix <- mutation_counts %>%
  pivot_wider(names_from = Variant_Classification, 
              values_from = mutation_count, 
              values_fill = 0) #%>%
  
# Step 10: Print the first few rows of the final matrix to verify that all columns are present
print(head(final_matrix_with_details))

# Step 11: Save the final matrix to a CSV file
write.csv(final_matrix_with_details, "/home/Ankita.Singh/Desktop/Computational framework identifies contributing factors for differences between PANoptosis clusters in multiple cancer types/neoantigent/complete_mutation_matrix_with_all_details.csv", row.names = FALSE)
