library(ggplot2)
library(dplyr)
library(ggpubr)
library(grid)

entropy_data <- read.csv('Smooth_std_entropy.csv', header = F)

colnames(entropy_data)[1] <- 'Genus'
entropy <- as.data.frame(t(entropy_data))
first_row <- entropy[1, ]
colnames(entropy) <- first_row
entropy <- entropy[-1,]
entropys <- as.data.frame(sapply(entropy, as.numeric))

base_position <- 1:1494
df1 <- cbind('Base_Position' = base_position, entropys)
df2 <- as.data.frame(df1)
df <- df2 %>% 
  dplyr::select(Base_Position, Actinoplanes, Azospirillum, Bacillus, Bradyrhizobium, Burkholderia, Ensifer,
                Enterobacter, Massilia, Mesorhizobium, Pseudomonas, Serratia, Streptomyces,
                Xanthomonas, Xylella)

# parameters
# region title offset
n = 0.018
# region title alpha
ra=0.8
# region title colour
rc='dark grey'
# panel alpha
pa=0.2
#panel colour
pc='grey'
#amplicon colour
ac = 'dark red'

genera <- colnames(df)[-1]

plots <- lapply(genera, function(genus) {
  
  plt <- ggplot(df)
  plt <- plt + geom_line(aes_string(x = "Base_Position", y = genus), lwd=1) +
    labs(title = genus) +
    theme(plot.title = element_text(face = "italic", size = 15))
  plt <- plt + theme_bw() + 
    scale_y_continuous(name = "", limits = c(-0.23, 1.06), 
                       breaks = c(0, 0.25, 0.5, 0.75, 1), 
                       labels = c("0.00", "0.25", "0.50", "0.75", "1.00"))
  
  # add colour panels
  plt <- plt + annotate('rect', xmin=6, xmax=104, ymin=0, ymax=1, fill=pc, alpha=pa)
  plt <- plt + annotate('rect', xmin=103, xmax=249, ymin=0, ymax=1, fill=pc, alpha=pa)
  plt <- plt + annotate('rect', xmin=331, xmax=501, ymin=0, ymax=1, fill=pc, alpha=pa)
  plt <- plt + annotate('rect', xmin=482, xmax=685, ymin=0, ymax=1, fill=pc, alpha=pa)
  plt <- plt + annotate('rect', xmin=749, xmax=879, ymin=0, ymax=1, fill=pc, alpha=pa)
  plt <- plt + annotate('rect', xmin=922, xmax=1031, ymin=0, ymax=1, fill=pc, alpha=pa)
  plt <- plt + annotate('rect', xmin=1060, xmax=1135, ymin=0, ymax=1, fill=pc, alpha=pa)
  plt <- plt + annotate('rect', xmin=1168, xmax=1326, ymin=0, ymax=1, fill=pc, alpha=pa)
  plt <- plt + annotate('rect', xmin=1331, xmax=1477, ymin=0, ymax=1, fill=pc, alpha=pa)
  
  # # add individual bars
  #V1
  plt <- plt + annotate('segment', x=6, xend=104, y=1, yend=1, colour=rc, alpha=ra)
  plt <- plt + annotate('text', x=6 + (104-9)/2, y=1+n+n, label='V1', cex=3)
  #V2
  plt <- plt + annotate('segment', x=103, xend=249, y=1+n, yend=1+n, colour=rc, alpha=ra)
  plt <- plt + annotate('text', x=103 + (249-103)/2, y=1+n+n+n, label='V2', cex=3)
  #V3
  plt <- plt + annotate('segment', x=331, xend=501, y=1, yend=1, colour=rc, alpha=ra)
  plt <- plt + annotate('text', x=331 + (501-331)/2, y=1+n+n, label='V3', cex=3)
  #V4
  plt <- plt + annotate('segment', x=482, xend=685, y=1+n, yend=1+n, colour=rc, alpha=ra)
  plt <- plt + annotate('text', x=482 + (685-482)/2, y=1+n+n+n, label='V4', cex=3)
  #V5
  plt <- plt + annotate('segment', x=749, xend=879, y=1, yend=1, colour=rc, alpha=ra)
  plt <- plt + annotate('text', x=749 + (879-749)/2, y=1+n+n, label='V5', cex=3)
  #V6
  plt <- plt + annotate('segment', x=922, xend=1031, y=1, yend=1, colour=rc, alpha=ra)
  plt <- plt + annotate('text', x=922 + (1031-922)/2, y=1+n+n, label='V6', cex=3)
  #V7
  plt <- plt + annotate('segment', x=1060, xend=1135, y=1, yend=1, colour=rc, alpha=ra)
  plt <- plt + annotate('text', x=1060 + (1135-1060)/2, y=1+n+n, label='V7', cex=3)
  #V8
  plt <- plt + annotate('segment', x=1168, xend=1326, y=1+n, yend=1+n, colour=rc, alpha=ra)
  plt <- plt + annotate('text', x=1168 + (1326-1168)/2, y=1+n+n+n, label='V8', cex=3)
  #V9
  plt <- plt + annotate('segment', x=1331, xend=1477, y=1, yend=1, colour=rc, alpha=ra)
  plt <- plt + annotate('text', x=1331 + (1477-1331)/2, y=1+n+n, label='V9', cex=3)
  
  # Add amplicon regions
  #V1-V3
  plt <- plt + annotate('segment', x=6, xend=501, y=-0.02-n, yend=-0.02-n, colour=ac, alpha=ra)
  plt <- plt + annotate('text', x=6 + (501-6)/2, y=-0.02-3*n, label='V1-V3', colour=ac, cex=3)
  #V3-V4
  plt <- plt + annotate('segment', x=331, xend=771, y=-0.03-2*n, yend=-0.03-2*n, colour=ac, alpha=ra)
  plt <- plt + annotate('text', x=331 + (771-331)/2, y=-0.03-3.5*n, label='V3-V4', colour=ac, cex=3)
  #V4-V5
  plt <- plt + annotate('segment', x=482, xend=886, y=-0.05-8*n, yend=-0.05-8*n, colour=ac, alpha=ra)
  plt <- plt + annotate('text', x=482 + (886-482)/2, y=-0.05-9.5*n, label='V4-V5', colour=ac, cex=3)
  #V4
  plt <- plt + annotate('segment', x=482, xend=771, y=-0.04-5*n, yend=-0.04-5*n, colour=ac, alpha=ra)
  plt <- plt + annotate('text', x=482 + (771-482)/2, y=-0.04-6.5*n, label='V4', colour=ac, cex=3)
  #V6-V8
  plt <- plt + annotate('segment', x=877, xend=1362, y=-0.02-n, yend=-0.02-n, colour=ac, alpha=ra)
  plt <- plt + annotate('text', x=877 + (1362-877)/2, y=-0.02-3*n, label='V6-V8', colour=ac, cex=3)
  #V6-V9
  plt <- plt + annotate('segment', x=877, xend=1466, y=-0.04-5*n, yend=-0.04-5*n, colour=ac, alpha=ra)
  plt <- plt + annotate('text', x=877 + (1466-877)/2, y=-0.04-6.5*n, label='V6-V9', colour=ac, cex=3)
  
  plt <- plt + theme(panel.grid = element_blank(),
                     plot.title = element_text(face = "italic", size = 15), 
                     axis.title.x = element_blank(), 
                     axis.title.y = element_blank())
  return(plt)
})

figure <- ggarrange(plotlist = plots, labels = c("a)","b)","c)","d)","e)","f)","g)","h)","i)","j)","k)","l)","m)","n)"), ncol = 3, nrow = 5)
figure1 <- annotate_figure(figure, 
                           left = textGrob("Standardized entropy", rot = 90, vjust = 1, gp = gpar(cex = 2)),
                           bottom = textGrob("Position along 16S gene", gp = gpar(cex = 2)))

ggsave("Supp_fig1.png", plot = figure1, width = 20, height = 25)