library(cowplot)
library(ggplot2)
library(ggtree)
library(TreeDist)
library(pvclust)
library(readr)
library(ape)
library(tidyr)
library(dplyr)

my_m <- read_tsv("ANIclustermap_matrix.tsv")
result <- pvclust(my_m, method.dist="cor", method.hclust="average", nboot=1000, parallel=TRUE)
tree <- as.phylo(result$hclust)

T1 <- ggtree(result)+
  geom_text(data=result, aes(label=bp), size = 2, hjust = -0.3)+ #verzija 2 za dodajanje bootstrap values
  geom_tiplab(align = TRUE, size = 2) +
  scale_y_reverse()+
  coord_cartesian(clip = "off")
T1

myfile2 <- file.path("Cupriavidus_MSA_V1V3_bt.fasta.treefile")
tree2 <- ape::read.tree(myfile2)
tree2 <- ape::root(tree2, "82633.13")
#write.tree(tree2, file = "Cupriavidus_V1V3_bt.newick")
T2 <- ggtree(tree2) +
  geom_text2(aes(subset = !isTip, label=label), size = 2, hjust = -0.3) + 
  geom_tiplab(hjust =1, align = TRUE, size = 2)+
  scale_x_reverse()+
  scale_y_reverse()+
  coord_cartesian(clip = "off")
T2

JRF_Cupriavidus <- JaccardRobinsonFoulds(tree, tree2, normalize = TRUE)
JRF <- round((JRF_Cupriavidus),3)


MCI_Cupriavidus <- MutualClusteringInfo(tree, tree2, normalize = TRUE)
MCI <- round((MCI_Cupriavidus),3)

dist_label <- paste("JRF =", JRF, "MCI =", MCI)

d1 = T1$data[T1$data$isTip,]
d1 <- d1 %>% select(-si, -au, -bp)
d1$x[] = 1  
d2 = T2$data[T2$data$isTip,]  
d2$x[] = 2  

TTcon <- rbind(d1, d2)  

L1 = ggplot(TTcon, aes(x = x, y = y, colour = label, group = label)) + geom_line(size = 1.1, color = "#5C5C5C") +   
  theme_void() + theme(legend.position="none", plot.margin = unit(c(1,0,1,0),"cm"))  

ANI_V1V3 <- cowplot::plot_grid(T1, L1 ,T2, labels = c("ANI dendrogram",dist_label, "V1-V3 region tree"),label_y = 0.93, nrow = 1, align = "hv", label_size = 10, scale = 0.9)
ANI_V1V3

ggsave("CupriavidusV1V3_1d.png", ANI_V1V3, width = 10, height = 6)

#V6V8 tree

myfile2 <- file.path("Cupriavidus_MSA_V6V8_bt.fasta.treefile")
tree2 <- ape::read.tree(myfile2)
tree2 <- ape::root(tree2, "82633.13")
T2 <- ggtree(tree2) +
  geom_text2(aes(subset = !isTip, label=label), size = 2, hjust = -0.3) + 
  geom_tiplab(hjust =1, align = TRUE, size = 2)+
  scale_x_reverse()+
  scale_y_reverse()+
  coord_cartesian(clip = "off")
T2

JRF_Cupriavidus <- JaccardRobinsonFoulds(tree, tree2, normalize = TRUE)
JRF <- round((JRF_Cupriavidus),3)


MCI_Cupriavidus <- MutualClusteringInfo(tree, tree2, normalize = TRUE)
MCI <- round((MCI_Cupriavidus),3)

dist_label <- paste("JRF =", JRF, "MCI =", MCI)

d1 = T1$data[T1$data$isTip,]
d1 <- d1 %>% select(-si, -au, -bp)
d1$x[] = 1  
d2 = T2$data[T2$data$isTip,]  
d2$x[] = 2  

TTcon <- rbind(d1, d2)  

L1 = ggplot(TTcon, aes(x = x, y = y, colour = label, group = label)) + geom_line(size = 1.1, color = "#5C5C5C") +   
  theme_void() + theme(legend.position="none", plot.margin = unit(c(1,0,1,0),"cm"))  

ANI_V6V8 <- cowplot::plot_grid(T1, L1 ,T2, labels = c("ANI dendrogram",dist_label, "V6-V8 region tree"),label_y = 0.93, nrow = 1, align = "hv", label_size = 10, scale = 0.9)
ANI_V6V8

ggsave("CupriavidusV6V8_1d.png", ANI_V6V8, width = 10, height = 6)

