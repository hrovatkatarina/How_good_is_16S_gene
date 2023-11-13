library(ggplot2)
library(tidyr)
library(dplyr)
theme_set(theme_bw())

distance_tableJRF <- read.csv('JRF_distance_ANI.csv')
distance_tableJRF_ANI10 <- read.csv('Distance_ANI_groups.csv')

JRF_selected <- distance_tableJRF %>% 
  filter(Genus == "Bacillus" | Genus == "Bradyrhizobium" | Genus == "Pseudomonas" | Genus == "Rhizobium" | 
           Genus == "Streptomyces") %>% 
  filter(Region == "Full_16S") %>% 
  mutate(Distancee = c("All genomes","All genomes","All genomes","All genomes","All genomes"))

JRF_selectedANI <- distance_tableJRF_ANI10 %>% 
  mutate(Distancee = c("One genome per ANI group"))

df_strain <- rbind(JRF_selected, JRF_selectedANI)

df_boxplot_ANI <- df_strain %>% 
  ggplot(aes(x = Distancee, y = Distance)) +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(aes(colour = Genus), width = 0.2, height = 0, alpha = 0.8, size = 2.5) +
  scale_color_manual(values = c("Bacillus" = "#458B74","Bradyrhizobium" = "#EE799F", "Pseudomonas" = "#1874CD", "Rhizobium" = "#984ea3", "Streptomyces" = "#FFD700")) +
  #geom_text(data = labels_wil, aes(x = Region, y = pos, label = Letters_wilcox), vjust = -0.1) +
  ylab("JRF distance: full 16S tree to ANI dendrogram")+
  xlab("")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 10),
        legend.text=element_text(size = 11),
        legend.title = element_text(size = 13),
        legend.position = c(0.85, 0.82))
ggsave("Boxplot_1f.png", df_boxplot_ANI, width = 6.8, height = 6)