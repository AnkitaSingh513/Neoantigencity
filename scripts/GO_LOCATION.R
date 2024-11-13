# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

# Load biomaRt library
library(biomaRt)

# Use Ensembl BioMart to fetch subcellular location data
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Load the list of genes from the CSV file
gene_data <- read.csv("Results/breast_cancer_gene_mutation_info.csv")
gene_list <- gene_data$Hugo_Symbol  # Extract the 'Hugo_Symbol' column

# Query BioMart to retrieve GO terms related to cellular components and their evidence codes
gene_data_biomart <- getBM(
  attributes = c('hgnc_symbol', 'go_id', 'name_1006', 'namespace_1003', 'go_linkage_type'),  # Use go_linkage_type for evidence codes
  filters = 'hgnc_symbol',
  values = gene_list,
  mart = ensembl
)

# Filter to include only subcellular location GO terms (cellular components)
cellular_component_data <- gene_data_biomart[gene_data_biomart$namespace_1003 == "cellular_component", ]

# Filter for high-confidence evidence codes (e.g., IDA, IMP, IPI)
high_confidence_codes <- c("IDA", "IMP", "IPI", "EXP", "TAS", "IC", "IEA")
accurate_go_terms <- cellular_component_data[cellular_component_data$go_linkage_type %in% high_confidence_codes, ]

# For merging: get only the top 2 GO terms per gene
cellular_component_aggregated <- aggregate(
  paste(go_id, name_1006, sep = ": ") ~ hgnc_symbol, 
  data = accurate_go_terms, 
  function(x) paste(head(x, 2), collapse = " | ")  # Limit to top 2 GO terms
)

# Rename the newly created column to be more descriptive
names(cellular_component_aggregated) <- c("Hugo_Symbol", "GO_Cellular_Component")

# Merge with the original gene data based on the 'Hugo_Symbol' column
merged_data <- merge(gene_data, cellular_component_aggregated, by = "Hugo_Symbol", all.x = TRUE)

# Save the updated file with accurate GO and subcellular location added at the end of the file
write.table(merged_data, "Results/breast_cancer_gene_summary_updated_with_GO_and_locations.csv", row.names = FALSE, quote=F, sep="\t")

# Print a preview of the updated data
print(head(merged_data))
