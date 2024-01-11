library(pvclust)
library(ape)
library(readr)

# This script will calculate bootstrap supported dendrogram of ANI matrix using pvclust. 


BT_dendro <- function(ANI_dir,ANI_file) {
  my_m <- read_tsv(ANI_file)
  result <- pvclust(my_m, method.dist="cor", method.hclust="average", nboot=1000, parallel=TRUE)
  phylo_tree <- as.phylo(result$hclust)
  bootstrap_values <- rev(result$edges$bp)*100
  phylo_tree$node.label <- paste0(phylo_tree$node.label, " (", bootstrap_values, ")")
  # Write the phylo object to a file in Newick format
  ANI_tree <- file.path(paste0(ANI_dir,"/ANI_tree_BT.newick"))
  write.tree(phylo_tree, file = ANI_tree)
 }
  
genera <- c("Actinoplanes", "Azospirillum", "Bacillus", "Bradyrhizobium", "Burkholderia", "Cupriavidus", "Ensifer",
            "Enterobacter", "Massilia", "Mesorhizobium", "Pseudomonas", "Rhizobium", "Serratia", "Streptomyces",
            "Xanthomonas", "Xylella")

for (genus in genera) {
  # Directory of ANI matrix
  ANI_dir <- file.path("Data/", genus)
  # ANI distance matrix
  ANI_file <- file.path(paste0(ANI_dir,"/ANIclustermap_result_new/ANIclustermap_matrix.tsv"))
  print(ANI_file)
  BT_dendro(ANI_dir, ANI_file)
}


