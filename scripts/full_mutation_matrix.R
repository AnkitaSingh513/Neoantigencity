# # Install necessary packages if not already installed
# if (!requireNamespace("maftools", quietly = TRUE)) {
#   install.packages("BiocManager")
#   BiocManager::install("maftools")
# }
# 
# if (!requireNamespace("tidyr", quietly = TRUE)) {
#   install.packages("tidyr")
# }

# Load required libraries
library(maftools)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
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

#Set the required path
setwd("/home/Raghvendra.Mall/TII/Projects/Raghav/CAR-T/Neoantigenicity/Neoantigencity/")

# Step 1: Load the MAF file (mutation data)
maf_file <- "/home/Raghvendra.Mall/TII/Projects/Raghav/CAR-T/General_Data/Ankita/TCGA_Data/mc3.v0.2.8.PUBLIC.maf.gz"
maf_data <- read.maf(maf = maf_file)

##################################################
# Step 2: Load the GDC sample file
gdc_file <- "Data/GDC-PANCAN.basic.phenotype.tsv"
gdc_data <- read_tsv(gdc_file)  # Use read_tsv for tab-delimited files

##################################################
# Step 3: Truncate 'sample' and 'Tumor_Sample_Barcode' to the first 12 characters
gdc_data <- gdc_data %>%
  mutate(sample_trimmed = substr(sample, 1, 12))  # Extract the first 12 characters from the GDC data

# Take GDC sample ids and match with maf data sample ids
gdc_data$sample_trimmed <- gdc_data$sample
maf_data@data$Tumor_Sample_Barcode <- as.character(as.vector(maf_data@data$Tumor_Sample_Barcode))
temp <- strsplit(maf_data@data$Tumor_Sample_Barcode,"-")
tumor_samples <- paste(unlist(map(temp,c(1))),unlist(map(temp,c(2))),unlist(map(temp,c(3))),unlist(map(temp,c(4))),sep = "-")

#Remove last character from samples
maf_data@data$Tumor_Short_Barcode <- str_sub(tumor_samples, end=-2)

# From MAF data get sample id, gene symbol, variant information, SNP id and chromosome no
mutation_data <- maf_data@data %>%
  select(Tumor_Short_Barcode, Hugo_Symbol, Variant_Classification, dbSNP_RS, Chromosome) %>%
  #mutate(Tumor_Sample_Barcode_trimmed = substr(Tumor_Sample_Barcode, 1, 12)) %>%
  filter(!is.na(Tumor_Short_Barcode), !is.na(Hugo_Symbol), !is.na(Variant_Classification))

######################################################
# Step 4: Merge the mutation data with GDC sample data
subset_gdc_data <- gdc_data[,c("sample_trimmed","cancer type abbreviation")]
colnames(subset_gdc_data) <- c("sample_trimmed","project_id")

merged_data <- mutation_data %>%
  left_join(subset_gdc_data, by = c("Tumor_Short_Barcode" = "sample_trimmed"))  %>% #, relationship = "many-to-many")
  filter(!is.na(project_id))

######################################################
# Step 5: Calculate global mutation frequencies
global_freq <- merged_data %>%
  group_by(Variant_Classification) %>%
  summarise(global_count = n()) #%>%
  #mutate(global_weight = global_count / sum(global_count))  # Normalize to get global weight

######################################################
## Step 6: Group by necessary columns and calculate mutation counts
mutation_counts <- merged_data %>%
  group_by(Tumor_Short_Barcode, Hugo_Symbol, Variant_Classification, dbSNP_RS, Chromosome, project_id) %>%
  summarise(mutation_count = n(), .groups = 'drop')

######################################################
# Step 7: Create a wide-format matrix including all specified mutation types as columns per sample per cancer
# Include all the required columns and mutation types
final_matrix <- mutation_counts %>%
  pivot_wider(names_from = Variant_Classification, 
              values_from = mutation_count, 
              values_fill = 0) 

#####################################################
# Step 8: Save the final matrix to a CSV file
write.csv(final_matrix, "Results/pancancer_per_sample_mutation_information.csv", row.names = FALSE, quote=F)

#####################################################
#Step 9: Get global summary of mutations for each gene
global_gene_summary <- maf_data@gene.summary
rev_global_gene_summary <- global_gene_summary %>%
  rowwise() %>%
  mutate(log_difference = calculate_difference_and_log(total, AlteredSamples))

rev_global_gene_summary <- rev_global_gene_summary[order(rev_global_gene_summary$log_difference,decreasing = T),]
write.csv(rev_global_gene_summary, file="Results/pancancer_per_gene_mutation_information.csv", row.names=F, quote=F)
