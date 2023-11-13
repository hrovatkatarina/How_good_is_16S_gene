library(tibble)
library(pheatmap)
library(tidyr)
library(ggplot2)

distance_table <- read.csv('JRF_distance_ANI.csv')

distance_tab <- distance_table %>% 
  pivot_wider(names_from = Region, values_from = Distance) %>% 
  column_to_rownames(var="Genus")

df_zscore <- apply(distance_tab, 1, function(x) (x - mean(x)) / sd(x))
df_zscore <- t(df_zscore)

sorted_rows <- c("Bacillus", "Streptomyces", "Actinoplanes", "Azospirillum", "Bradyrhizobium", "Mesorhizobium", 
                 "Rhizobium", "Ensifer", "Massilia", "Cupriavidus", "Burkholderia","Xylella", "Xanthomonas", 
                 "Pseudomonas", "Serratia", "Enterobacter")

dfJRF <- as.matrix(df_zscore)
dfJRF_sorted <- dfJRF[match(sorted_rows, rownames(dfJRF)), , drop = FALSE]

Class <- data.frame(Class = c("Bacilli", "Actinomycetia", "Actinomycetia", "Alphaproteobacteria", 
                              "Alphaproteobacteria", "Alphaproteobacteria", "Alphaproteobacteria", 
                              "Alphaproteobacteria", "Gammaproteobacteria", "Gammaproteobacteria", 
                              "Gammaproteobacteria", "Gammaproteobacteria","Gammaproteobacteria", 
                              "Gammaproteobacteria", "Gammaproteobacteria", "Gammaproteobacteria"))

row.names(Class) <- row.names(dfJRF_sorted)

order_color <- list(Class = c(Bacilli = "#528B8B", Actinomycetia = "#e41a1c", Alphaproteobacteria ="#CD96CD", Gammaproteobacteria = "#A2CD5A"))

plt1JRF <- pheatmap(dfJRF_sorted,
                    main = "Jaccard-Robinson-Foulds distance: Full length 16S gene\nand variable regions tree vs. ANI dendrogram",
                    cluster_cols = F,
                    cluster_rows = T,
                    cellwidth = 30,
                    cellheight = 20,
                    border_color = NA,
                    gaps_row = 1:16,
                    color = rev(hcl.colors(50, "Reds")),
                    annotation_colors = order_color,
                    annotation_row = Class,
                    legend_breaks = c(1.87, 1, 0, -1, -2),
                    legend_labels = c("Row Z-score", 1, 0, -1, -2))
ggsave("Heatmap_1e.png", plt1JRF, width = 9, height = 7)