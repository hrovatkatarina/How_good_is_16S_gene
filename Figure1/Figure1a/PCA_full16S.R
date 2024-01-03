library(ggplot2)
library(ggfortify)
library(Biostrings)

theme_set(theme_bw())

fasta_file <- "full16S_MSAall.fasta"
fasta_text <- readLines(fasta_file)
fixed_fasta_text <- gsub("U", "T", fasta_text)
fixed_fasta_file <- tempfile() 
writeLines(fixed_fasta_text, fixed_fasta_file)


# Import multiple sequence alignment of full-length 16S rRNA gene 
dna_seqs <- readDNAStringSet(fixed_fasta_file)

# Calculate oligonulcleotide frequency
kmer_freq_nonames <- oligonucleotideFrequency(dna_seqs, width=4, step=1,
                                              as.prob=FALSE, as.array=FALSE,
                                              fast.moving.side="right", with.labels=TRUE)

# Reassign names to the resulting kmer_freq_nonames
Genus <- sapply(strsplit(names(dna_seqs)," "), `[`, 2)
name_df <- as.data.frame(Genus)

# PCA
kmer.pca <- prcomp(kmer_freq_nonames,
                   center = TRUE,
                   scale. = TRUE)

kmer.pca.plot <- autoplot(kmer.pca, data = name_df, colour = 'black', shape = 'Genus', fill = "Genus", size = 2.5) +
  ggtitle("Oligonucleotide frequencies across full 16S rRNA gene")+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.text=element_text(size = 11, hjust = 0, face="italic"),
        plot.title = element_text(size = 14.5))

kmer.pca.plot_final <- kmer.pca.plot + 
  scale_fill_manual(values = c("Bacillus" = "#528B8B","Streptomyces"="#FFD700","Actinoplanes" = "#e41a1c",
                               "Azospirillum" = "#ff7f00", "Bradyrhizobium" = "#FFB5C5","Mesorhizobium" = "#CD96CD",
                               "Rhizobium" = "#984ea3","Ensifer" = "#FFBBFF", "Massilia" = "#4daf4a","Cupriavidus" = "#A2CD5A",
                               "Burkholderia" = "#B4EEB4","Pseudomonas" = "#377eb8","Serratia" = "#a65628",
                               "Enterobacter" = "#CDAA7D", "Xylella" = "#C7C7C7","Xanthomonas" = "#595959", 
                               "Paenibacillus" = "#528B8B","Kitasatospora" = "#FFD700","Micromonospora" = "#e41a1c",
                               "Nitrospirillum" = "#ff7f00","Rhodopseudomonas" = "#FFB5C5","Aminobacter" = "#CD96CD",
                               "Agrobacterium" = "#984ea3","Sinorhizobium" = "#FFBBFF","Duganella" = "#4daf4a", 
                               "Ralstonia" = "#A2CD5A","Paraburkholderia" = "#B4EEB4","Azotobacter" = "#377eb8",
                               "Rahnella" = "#a65628", "Citrobacter" = "#CDAA7D","Pseudoxanthomonas" = "#595959"))+
  scale_shape_manual(values = c("Bacillus" = 21, "Paenibacillus" = 24, "Streptomyces"=21, "Kitasatospora" = 24, 
                                "Actinoplanes" = 21, "Micromonospora" = 24, "Rhizobium" = 21,"Agrobacterium" = 24, 
                                "Mesorhizobium" = 21, "Aminobacter" = 24, "Ensifer" = 21, "Sinorhizobium" = 24, 
                                "Bradyrhizobium" = 21, "Rhodopseudomonas" = 24, "Massilia" = 21, "Duganella" = 24,
                                "Cupriavidus" = 21, "Ralstonia" = 24, "Burkholderia" = 21, "Paraburkholderia" = 24, 
                                "Enterobacter" = 21, "Citrobacter" = 24, "Serratia" = 21, "Rahnella" = 24, 
                                "Pseudomonas" = 21, "Azotobacter" = 24, "Xylella" = 21, "Xanthomonas" = 21, 
                                "Pseudoxanthomonas" =24, "Azospirillum" = 21, "Nitrospirillum" = 24))

kmer.pca.plot_final
ggsave(filename = "PCA_full16S.png", plot=kmer.pca.plot_final, width = 9, height = 5)