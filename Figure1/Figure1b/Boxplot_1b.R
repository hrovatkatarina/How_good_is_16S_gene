library(ggplot2)
library(ggpubr)
theme_set(theme_bw())

df_boxplot <- read.csv("Coregenome_distance.csv")

#Boxplot
df_box_plot <- df_boxplot %>% 
  ggplot(aes(x = Region, y = Distance)) +
  geom_boxplot(width = 0.4, outlier.shape = NA) +
  geom_jitter(aes(colour = Genus), width = 0.2, height = 0, alpha = 0.8, size = 2.5) +
  scale_color_manual(
     values = c("Actinoplanes" = "#e41a1c","Azospirillum" = "#ff7f00", "Cupriavidus" = "#a6d96a", "Ensifer" = "#984ea3", "Massilia" = "#1a9641"),
     labels = c(expression(italic("Actinoplanes")), expression(italic("Azospirillum")), expression(italic("Cupriavidus")), expression(italic("Ensifer")),
               expression(italic("Massilia")))) +
  ylab("JRF distance to single-copy marker genes tree")+
  xlab("")+
  coord_cartesian(ylim = c(0.1, 0.8)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(size = 13),
        axis.text.x = element_text(size = 10),
        legend.text=element_text(size = 11, hjust = 0),
        legend.title = element_text(size = 13),
        legend.position = c(0.87, 0.21)) +
  stat_compare_means(label = "p.signif", size = 7, label.y = 0.75, ref.group = "ANI")
ggsave("Boxplot_1b.png", df_box_plot, width = 6.8, height = 5)
