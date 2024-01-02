library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
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

df_strain2 <- df_strain %>% 
  filter(Distancee == "One genome per ANI group") %>%
  mutate(Shape = "One genome per ANI group")

level_order <- c("Pseudomonas", "Streptomyces", "Bacillus", "Rhizobium", "Bradyrhizobium")

df_dots <- df_strain %>% 
  filter(Distancee == "All genomes") %>% 
  mutate(Color = "All genomes")

df_boxplot_ANI2 <- df_strain2 %>% 
  ggplot(aes(x = factor(Genus, level = level_order), y = Distance, color = "All")) +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.8, size = 3.5, shape = 16) +
  geom_point(data = df_dots, aes(x = Genus, y = Distance, color = 'All genomes'), size = 4, shape = 17) +
  scale_color_manual(name='The set of genomes used to\ncalculate the phylogenetic tree:',
                     labels=c('One genome per ANI group', 'All genomes'),
                     values=c("black", "#8C8C8C")) +
  ylab("JRF distance: full 16S tree to ANI dendrogram")+
  xlab("")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 13, face = "italic"),
        legend.position = c(0.79, 0.88),
        legend.title = element_text(size = 12.5),
        legend.text=element_text(size=11.5)) +
  stat_compare_means(method = 'anova', label.x.npc = 0.01)
ggsave("Boxplot_1f.png", df_boxplot_ANI2, width = 6.8, height = 7)
df_boxplot_ANI2