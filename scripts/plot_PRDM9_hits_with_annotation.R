args <- commandArgs()
path_fimo_window_bed <- args[6]
x_min <- as.numeric(args[7])
x_max <- as.numeric(args[8])
num_chr <- as.numeric(args[9])
width <- as.numeric(args[10])
path_annotation <- args[11]
path_output <- args[12]

library(ggridges)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(scales) # for pretty_breaks()

options(scipen = 9)

x <- read.delim(path_fimo_window_bed, header = F)
colnames(x) <- c('chrom', 'ref.begin', 'ref.end', 'hits')

#x$chrom <- gsub('[:].*$','', x$chrom)

if (num_chr == 13) {
  colors <- c("#F8766D")
} else if (num_chr == 14) {
  colors <- c("#A3A500")
} else if (num_chr == 15) {
  colors <- c("#00BF7D")
} else if (num_chr == 21) {
  colors <- c("#00B0F6")
} else if (num_chr == 22) {
  colors <- c("#E76BF4")
}else {
  colors <- c("#000000")
}

chr <- paste0('chm13#chr', num_chr)
xx <- x[x$chrom == chr & x$ref.begin >= x_min & x$ref.end <= x_max,]
p <- ggplot(xx, aes(
  x = (ref.begin + (ref.end - ref.begin) / 2) / 1000000, width = ref.end - ref.begin,
  y=hits,
  color=chrom)
) +
  geom_step(size=0.3) +
  #facet_wrap(~chrom, scales = "free_y", ncol=1, labeller = labeller(variable = labels)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    
    text = element_text(size = 10),
    axis.text.x = element_text(size = 10, angle = 45, vjust = 1.0, hjust=1),
    axis.text.y = element_text(size = 10),
    
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.position = "none",
    
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    plot.margin = unit(c(0.1,0.1,0.1,1.025), "cm")
  ) + scale_x_continuous(
    #limits = c(x_min, x_max),
    breaks = pretty_breaks(n=20),
    expand = c(0.01, 0.01)) + labs(
      x = paste('Position (Mbp)'),
      y = paste('PRDM9 motif hits'),
      color = "Chromosome"
    ) + scale_y_continuous(
      #limits = c(0, max(10, max(xx$hits))),
      breaks=pretty_breaks(n=6)
    ) +
  #ggtitle('PRDM9 density hits for each bp in the SST1 repetitive unit') +
  scale_color_manual(values=colors) +
  guides(colour = guide_legend(override.aes = list(size=10)))
#p
#ggsave(plot = p, path_output, width = width, height = length(unique(x$chrom))*4, units = "cm", dpi = 300, bg = "transparent", limitsize = FALSE)


library(png)
library(grid)
img <- readPNG(path_annotation)

ggplotted_img <- ggplot() +
  annotation_custom(
    rasterGrob(img, width = 1, height = 1),
    xmin = - Inf, xmax = Inf,
    ymin = - Inf, ymax = Inf
  ) + theme(
    
    plot.margin = unit(c(0,0.4,0,-0.7), "cm")
  )

library(ggpubr)
p_with_annotation <- ggpubr::ggarrange(
  ggplotted_img, p,
  labels=c('', ''),
  heights = c(1.95, 3),
  legend = "none", # legend position,
  common.legend = T,
  nrow = 2
)

ggsave(
  plot = p_with_annotation,
  path_output,
  width = width, height = 8,
  units = "cm",
  dpi = 300, bg = "transparent",
  limitsize = FALSE
)
