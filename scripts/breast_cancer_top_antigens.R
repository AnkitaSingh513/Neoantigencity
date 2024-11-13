library(data.table)
library(ggplot2)
library(data.table)
library(tidyverse)
library(ggpicrust2)
library(dada2)
library(stringi)
library(ComplexHeatmap)
library(rstatix)
library(ggpubr)

setwd("/home/Raghvendra.Mall/TII/Projects/Raghav/CAR-T/Neoantigenicity/Neoantigencity/")

#Get dataset for breast neoantigenicity 
breast_neoantigen_df <- fread("Results/breast_cancer_gene_summary_updated_with_GO_and_locations.csv",header=T,sep="\t")
temp_annotations <- str_replace_na(unlist(lapply(strsplit(breast_neoantigen_df$GO_Cellular_Component, split=" & "),'[[',1)), replacement = "GO: X")
cellular_localization_info <- unlist(lapply(strsplit(temp_annotations,split=": "),'[[',2))

#Get the mutation matrix
order_ids <- order(breast_neoantigen_df$log_difference, decreasing=T)

breast_neoantigen_df <- breast_neoantigen_df[order_ids,]
aggr_mut_matrix <- breast_neoantigen_df %>% select(altered, total, Missense_Mutation, Frame_Shift_Del, Nonsense_Mutation, Splice_Site, In_Frame_Del, Frame_Shift_Ins, Translation_Start_Site, Nonstop_Mutation, In_Frame_Ins) %>% as.matrix()
rownames(aggr_mut_matrix) <- breast_neoantigen_df$Hugo_Symbol

#Get the ranking score
aggr_ranking_score <- breast_neoantigen_df$log_difference

subset_aggr_mut_matrix <- aggr_mut_matrix[c(1:50),]
subset_aggr_ranking_score <- aggr_ranking_score[c(1:50)]
rev_cellular_localization_info <- cellular_localization_info[order_ids]
subset_cellular_localization_info <- rev_cellular_localization_info[c(1:50)]
subset_cellular_localization_info[subset_cellular_localization_info=="cytoplasmic ribonucleoprotein granule"] <- "cytoplasm"

#Make the complex Heatmap
library(RColorBrewer)
library(circlize)
rank_col_fun = colorRamp2(c(3,8), c("blue","red")) 
local_col_fun = RColorBrewer::brewer.pal(n=10,name="Set3")
col_fun = colorRamp2(breaks=c(0,1,10,100,1000), colors = c("blue","cyan","white","maroon","red"))

row_ha = rowAnnotation(Rank = subset_aggr_ranking_score,
                       Localization = subset_cellular_localization_info,
                       show_legend = T,
                       col = list(Rank = rank_col_fun,
                                  Localization = c("anchoring junction"=local_col_fun[1],
                                                   "cell projection"=local_col_fun[2],
                                                   "cytoplasm"=local_col_fun[3],
                                                   "cytosol"=local_col_fun[4],
                                                   "dynein complex"=local_col_fun[5],
                                                   "extracellular region"=local_col_fun[6],
                                                   "membrane"=local_col_fun[7],
                                                   "nuclear outer membrane"=local_col_fun[8],
                                                   "nucleoplasm"=local_col_fun[9],
                                                   "nucleus"=local_col_fun[10],
                                                   "X"=local_col_fun[11])))

ht <- Heatmap(subset_aggr_mut_matrix, 
              name="Mutations / Samples",
              cluster_rows=F,
              col = col_fun,
              width = unit(5, "in"), 
              height = unit(10, "in"),
              rect_gp = gpar(col = "white", lwd = 2),
              column_dend_height = unit(1, "cm"),
              clustering_method_columns = "ward.D2", 
              row_names_gp = gpar(fontsize = 11, box_col = rep(c("red"), times = 50)),
              row_title_rot = 0,
              column_names_gp = gpar(fontsize=11),
              column_names_side = "bottom",
              column_title = NULL,
              right_annotation = row_ha,
)

pdf("Results/Breast_Cancer_Top_Antigens_with_Annotations_Heatmap.pdf", height = 10, width=6)
draw(ht)
dev.off()

