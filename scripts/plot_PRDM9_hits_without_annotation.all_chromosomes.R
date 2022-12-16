args <- commandArgs()
path_fimo_window_bed <- args[6]
scale <- as.numeric(args[7])
x_axis_label <- args[8]
chrom_suffix <- args[9]
width <- as.numeric(args[10])
label_color <- args[11]
path_output <- args[12]

library(ggridges)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(scales) # for pretty_breaks()

options(scipen = 9)

x <- read.delim(path_fimo_window_bed, header = F)
colnames(x) <- c('chrom', 'ref.begin', 'ref.end', 'hits')

x$group <- gsub('[:].*$', chrom_suffix, x$chrom)
xx <- x

if(length(unique(xx$chrom)) > 3) {
  colors <- c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF4")
} else {
  # Assume 3 (chr13/chr14/chr21)
  colors <- c("#F8766D", "#A3A500", "#00B0F6")
}

p <- ggplot(xx, aes(
  x = (ref.begin + (ref.end - ref.begin) / 2) / scale, width = ref.end - ref.begin,
  y=hits,
  color=group)
) +
  geom_step() +
  facet_wrap(~chrom, scales = "free_y", ncol=1)+#, labeller = labeller(variable = labels)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    
    text = element_text(size = 18),
    axis.text.x = element_text(size = 12, angle = 45, vjust = 1.0, hjust=1),
    axis.text.y = element_text(size = 12),
    
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.position = "right",
    
    #strip.text.x = element_blank(),
    #strip.text.y = element_blank(),
    plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
  ) + scale_x_continuous(
    #limits = c(x_min, x_max),
    breaks = pretty_breaks(n=20),
    expand = c(0.01, 0.01)) + labs(
      x = paste(x_axis_label),
      y = paste('PRDM9 motif hits'),
      color = label_color
    ) + scale_y_continuous(
      #limits = c(0, max(10, max(xx$hits))),
      breaks=pretty_breaks(n=6)
    ) +
  #ggtitle('PRDM9 density hits for each bp in the SST1 repetitive unit') +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(override.aes = list(size=10)))
#p
height <- max(8, length(unique(x$chrom))*4.3)
ggsave(plot = p, path_output, width = width, height = height, units = "cm", dpi = 300, bg = "transparent", limitsize = FALSE)
