library(TreeDist)
library(ggtree)

# Distances: ANI, 16S and V-regions vs core-genome tree (for boxplot).

genera <- c("Actinoplanes", "Azospirillum", "Cupriavidus", "Ensifer", "Massilia")
distance_table <- data.frame(Genus = character(), Distance = numeric())
regions <- c("ANI", "Full 16S", "V1V3", "V3V4", "V4", "V4V5", "V6V8", "V6V9")

for (genus in genera) {
  my_dir <- file.path("C:/Users/Uporabnik/Desktop/Wageningen/Tree", genus)
  
  myfile <- file.path(my_dir, "core_gene_alignment.aln.treefile")
  tree_ref <- ape::read.tree(myfile)
  
  myfile <- file.path(my_dir, "ANI_tree_BT.newick")
  tree_ANI <- ape::read.tree(myfile)
  
  myfile <- file.path(my_dir, paste0(genus,"_full16s.fasta.treefile"))
  tree_full16S <- ape::read.tree(myfile)

  myfile <- file.path(my_dir, paste0(genus,"_MSA_V1V3.fasta.treefile"))
  tree_1 <- ape::read.tree(myfile)

  myfile <- file.path(my_dir, paste0(genus,"_MSA_V3V4.fasta.treefile"))
  tree_3 <- ape::read.tree(myfile)

  myfile <- file.path(my_dir, paste0(genus,"_MSA_V4.fasta.treefile"))
  tree_4 <- ape::read.tree(myfile)

  myfile <- file.path(my_dir, paste0(genus,"_MSA_V4V5.fasta.treefile"))
  tree_5 <- ape::read.tree(myfile)

  myfile <- file.path(my_dir, paste0(genus,"_MSA_V6V8.fasta.treefile"))
  tree_6 <- ape::read.tree(myfile)

  myfile <- file.path(my_dir, paste0(genus,"_MSA_V6V9.fasta.treefile"))
  tree_9 <- ape::read.tree(myfile)
  
  trees <- list(tree_ANI, tree_full16S, tree_1, tree_3, tree_4, tree_5, tree_6, tree_9)
  
  distance <- JaccardRobinsonFoulds(tree_ref, trees, normalize = TRUE)
  
  distance_table <- rbind(distance_table, data.frame(Genus = genus, Distance = distance))
}

distance_table <- cbind(distance_table, Region = regions)


write.csv(distance_table, "C:\\Users\\Uporabnik\\Desktop\\Wageningen\\Tree\\Coregenome_distance.csv")

# Distances: Full 16S, V-regions for all genera vs. ANI (for heatmap)

genera <- c("Actinoplanes", "Azospirillum", "Bacillus", "Bradyrhizobium", "Burkholderia", "Cupriavidus", "Ensifer", 
            "Enterobacter", "Massilia", "Mesorhizobium", "Pseudomonas", "Rhizobium", "Serratia", "Streptomyces", 
            "Xanthomonas", "Xylella")

distance_table_outputJRF <- data.frame(Genus = character(), Distance = numeric())
distance_tableJRF <- data.frame(Genus = character(), Distance = numeric())

#distance_table_outputMCI <- data.frame(Genus = character(), Distance = numeric())
#distance_tableMCI <- data.frame(Genus = character(), Distance = numeric())

regions <- c("Full 16S", "V1V3", "V3V4", "V4", "V4V5", "V6V8", "V6V9")

for (genus in genera) {
  my_dir <- file.path("C:\\Users\\Uporabnik\\Desktop\\Wageningen\\Tree", genus)
  myfile <- file.path(paste0(my_dir,"\\ANI_tree_BT.newick"))
  tree_ANI <- ape::read.tree(myfile)
  
  region = "_full16S"
  myfile <- file.path(my_dir, paste0(genus, region, ".fasta.treefile"))
  tree_16s <- ape::read.tree(myfile)
  
  region1 <- "_MSA_V1V3"
  myfile <- file.path(my_dir, paste0(genus, region1, ".fasta.treefile"))
  tree_1 <- ape::read.tree(myfile)
  
  region3 <- "_MSA_V3V4"
  myfile <- file.path(my_dir, paste0(genus, region3, ".fasta.treefile"))
  tree_3 <- ape::read.tree(myfile)
  
  region4 <- "_MSA_V4"
  myfile <- file.path(my_dir, paste0(genus, region4, ".fasta.treefile"))
  tree_4 <- ape::read.tree(myfile)
  
  region5 <- "_MSA_V4V5"
  myfile <- file.path(my_dir, paste0(genus, region5, ".fasta.treefile"))
  tree_5 <- ape::read.tree(myfile)
  
  region6 <- "_MSA_V6V8"
  myfile <- file.path(my_dir, paste0(genus, region6, ".fasta.treefile"))
  tree_6 <- ape::read.tree(myfile)
  
  region9 <- "_MSA_V6V9"
  myfile <- file.path(my_dir, paste0(genus, region9, ".fasta.treefile"))
  tree_9 <- ape::read.tree(myfile)
  
  trees <- list(tree_16s, tree_1, tree_3, tree_4, tree_5, tree_6, tree_9)
  
  JRF <- JaccardRobinsonFoulds(tree_ANI, trees, normalize = TRUE)
  #MCI <- MutualClusteringInfo(tree_ANI, trees, normalize = TRUE)
  
  distance_table_outputJRF <- rbind(distance_table_outputJRF, data.frame(Genus = genus, Distance = JRF))
  #MCI_output <- rbind(distance_table_outputMCI, data.frame(Genus = genus, Distance = MCI))
  
}

distance_tableJRF <- cbind(distance_table_outputJRF, Region = regions)

#tableMCI <- cbind(MCI_output, Region = regions)

write.csv(distance_tableJRF, 'JRF_disatnce_ANI.csv', row.names = FALSE)
#write.csv(tableMCI, 'MCI_ANI.csv')

# Distances: 16S and one per ANI group vs ANI (for genera with a lot of genomes)




