# Load required libraries
library(maftools)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# Define a function to calculate the difference and log for each gene
calculate_difference_and_log <- function(total, mutated) {
  # Subtract 'MutatedSamples' from 'total'
  if (is.na(total) | is.na(mutated)) {
    return(NA)
  }
  
  diff <- (total / mutated)* log(mutated+1)
  
  return(diff)
}

# Set working directory (optional, adjust to your path)
setwd("/home/Raghvendra.Mall/TII/Projects/Raghav/CAR-T/Neoantigenicity/Neoantigencity/")

# Step 1: Load the MAF file (mutation data)
maf_file <- "/home/Raghvendra.Mall/TII/Projects/Raghav/CAR-T/General_Data/Ankita/TCGA_Data/mc3.v0.2.8.PUBLIC.maf.gz"
maf_data <- read.maf(maf = maf_file)

# Step 2: Load the GDC sample file
gdc_file <- "Data/GDC-PANCAN.basic.phenotype.tsv"
gdc_data <- read_tsv(gdc_file)  # Read GDC data

# Step 3: Clean and prepare 'sample' and 'Tumor_Sample_Barcode'
gdc_data <- gdc_data %>%
  mutate(sample_trimmed = substr(sample, 1, 12))  # Extract the first 12 characters for matching

# Ensure `gdc_data` has unique `sample_trimmed` entries to avoid duplicate rows in join
gdc_data <- gdc_data %>%
  distinct(sample_trimmed, .keep_all = TRUE)

# Create shortened tumor sample barcodes
maf_data@data$Tumor_Short_Barcode <- substr(maf_data@data$Tumor_Sample_Barcode, 1, 12)

# Extract relevant columns for mutation data
mutation_data <- maf_data@data %>%
  select(Tumor_Short_Barcode, Hugo_Symbol, Variant_Classification, dbSNP_RS, Chromosome) %>%
  filter(!is.na(Tumor_Short_Barcode), !is.na(Hugo_Symbol), !is.na(Variant_Classification))
head(mutation_data)

# Step 4: Merge the mutation data with GDC sample data based on the short barcode
subset_gdc_data <- gdc_data[,c("sample_trimmed","cancer type abbreviation")]
colnames(subset_gdc_data) <- c("sample_trimmed","project_id")

merged_data <- mutation_data %>%
  left_join(subset_gdc_data, by = c("Tumor_Short_Barcode" = "sample_trimmed"), relationship = "many-to-many") %>%
  filter(!is.na(project_id))  # Ensure only samples with a valid project_id are retained
head(merged_data)

# Step 5: Calculate mutation counts grouped by relevant columns
mutation_counts <- merged_data %>%
  group_by(Hugo_Symbol, Variant_Classification, Tumor_Short_Barcode, project_id) %>%
  summarise(mutation_count = n(), .groups = 'drop')

altered_samples_count_per_cancer_per_gen <- mutation_counts %>% group_by(Hugo_Symbol, project_id)  %>% summarise(altered = n(), .groups="drop") 

total_mutations_count_per_cancer_per_gen <- merged_data %>% group_by(Hugo_Symbol, project_id) %>% summarise(total = n(), .groups="drop") 

merged_counts <- left_join(altered_samples_count_per_cancer_per_gen,total_mutations_count_per_cancer_per_gen,by = c("Hugo_Symbol", "project_id"))  
 
head(merged_counts)

merged_counts <- left_join(altered_samples_count_per_cancer_per_gen, 
                           total_mutations_count_per_cancer_per_gen, 
                           by = c("Hugo_Symbol", "project_id"))

final_merged_data <- left_join(merged_counts, 
                               mutation_counts, 
                               by = c("Hugo_Symbol", "project_id"))
head(final_merged_data)


temp_samples_count_per_cancer_per_gen <-  final_merged_data %>% group_by(Hugo_Symbol, project_id,altered,total,Variant_Classification)  %>% summarise(mutation_count = n(), .groups="drop")

# Step 6: Create a wide-format matrix to include mutation types as columns
final_matrix <- temp_samples_count_per_cancer_per_gen %>%
  pivot_wider(names_from = Variant_Classification, 
              values_from = mutation_count, 
              values_fill = 0)

gene_summary_per_cancer <- final_matrix %>%
  rowwise() %>%
  mutate(log_difference = calculate_difference_and_log(total, altered))  # Use correct column names here

write.csv(gene_summary_per_cancer,"Results/per_cancer_gene_mutation_info.csv",row.names=F, quote=F)

#Do the analysis for breast cancer
################################################################################
breast_cancer_gene_summary <- gene_summary_per_cancer[gene_summary_per_cancer$project_id=="BRCA",] %>% data.frame()
breast_cancer_gene_summary <- breast_cancer_gene_summary[order(breast_cancer_gene_summary$log_difference,decreasing = T),]
head(breast_cancer_gene_summary, n = 20)
print(breast_cancer_gene_summary)

# Save the breast_cancer_gene_summary to a CSV file
write.csv(breast_cancer_gene_summary, "Results/breast_cancer_gene_mutation_info.csv", row.names = FALSE, quote=F)



# # Calculate `total` as the sum of all mutation types for each row
# final_matrix <- final_matrix %>%
#   mutate(total = rowSums(select(., -c(Hugo_Symbol, Tumor_Short_Barcode, project_id))))
# head(final_matrix)
# # Calculate the correct `MutatedSamples` and `AlteredSamples` per sample and gene
# # Group by both `Hugo_Symbol` and `Tumor_Short_Barcode` for accurate calculation
# # Correct the calculation for `MutatedSamples`
# final_matrix <- final_matrix %>%
#   mutate(
#     MutatedSamples = ifelse(total > 0, total, 0)  # MutatedSamples should match `total` if any mutations exist
#   )
# 
# # Calculate `AlteredSamples` for each `Hugo_Symbol`
# gene_level_summary <- final_matrix %>%
#   group_by(Hugo_Symbol) %>%
#   summarise(AlteredSamples = n_distinct(Tumor_Short_Barcode[total > 0]))  # Count unique samples with mutations
# 
# # Drop any existing `AlteredSamples` columns before merging to avoid conflicts
# final_matrix <- final_matrix %>% select(-starts_with("AlteredSamples"))
# 
# # Merge `AlteredSamples` back into `final_matrix`
# final_matrix <- final_matrix %>%
#   left_join(gene_level_summary, by = "Hugo_Symbol")
# 
# # Reorder columns to match expected output
# final_matrix <- final_matrix %>%
#   select(Hugo_Symbol, Tumor_Short_Barcode, project_id, total, MutatedSamples, AlteredSamples, everything())
# 
# # Print the first few rows of the final matrix to verify the output
# print(head(final_matrix))
# 
# # Save the final matrix to a CSV file
# output_file <- "/home/Ankita.Singh/Desktop/Computational framework identifies contributing factors for differences between PANoptosis clusters in multiple cancer types/neoantigent/complete_mutation_matrix_with_all_details.csv"
# write.csv(final_matrix, output_file, row.names = FALSE)