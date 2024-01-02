library(TreeDist)
library(ggtree)

# Distances: ANI, 16S and V-regions vs core-genome tree (for boxplot).

genera <- c("Actinoplanes", "Azospirillum", "Cupriavidus", "Ensifer", "Massilia")
distance_table <- data.frame(Genus = character(), Distance = numeric())
regions <- c("ANI", "Full 16S", "V1V3", "V3V4", "V4", "V4V5", "V6V8", "V6V9")

for (genus in genera) {
  
  myfile <- file.path(genus, "core_gene_alignment.aln.treefile")
  tree_ref <- ape::read.tree(myfile)
  
  myfile <- file.path(genus, "ANI_tree_BT.newick")
  tree_ANI <- ape::read.tree(myfile)
  
  myfile <- file.path(genus, paste0(genus,"_full16s.fasta.treefile"))
  tree_full16S <- ape::read.tree(myfile)

  myfile <- file.path(genus, paste0(genus,"_MSA_V1V3.fasta.treefile"))
  tree_1 <- ape::read.tree(myfile)

  myfile <- file.path(genus, paste0(genus,"_MSA_V3V4.fasta.treefile"))
  tree_3 <- ape::read.tree(myfile)

  myfile <- file.path(genus, paste0(genus,"_MSA_V4.fasta.treefile"))
  tree_4 <- ape::read.tree(myfile)

  myfile <- file.path(genus, paste0(genus,"_MSA_V4V5.fasta.treefile"))
  tree_5 <- ape::read.tree(myfile)

  myfile <- file.path(genus, paste0(genus,"_MSA_V6V8.fasta.treefile"))
  tree_6 <- ape::read.tree(myfile)

  myfile <- file.path(genus, paste0(genus,"_MSA_V6V9.fasta.treefile"))
  tree_9 <- ape::read.tree(myfile)
  
  trees <- list(tree_ANI, tree_full16S, tree_1, tree_3, tree_4, tree_5, tree_6, tree_9)
  
  distance <- JaccardRobinsonFoulds(tree_ref, trees, normalize = TRUE)
  
  distance_table <- rbind(distance_table, data.frame(Genus = genus, Distance = distance))
}

distance_table <- cbind(distance_table, Region = regions)


write.csv(distance_table, "Coregenome_distance.csv")



