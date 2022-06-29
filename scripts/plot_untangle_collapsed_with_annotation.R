args <- commandArgs()
path_support_dedup_tsv <- args[6]
x_min <- as.numeric(args[7])
x_max <- as.numeric(args[8])
width <- as.numeric(args[9])
num_chr <- as.numeric(args[10])
path_annotation <- args[11]
prefix_output <- args[12]


# library
library(ggridges)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(scales) # for pretty_breaks()

x <- read.delim(path_support_dedup_tsv)
colnames(x) <- c("ground.target", "start", "end", "chr13", "chr14", "chr15", "chr21", "chr22")

# Apply filters
chr <- paste0('chm13#chr', num_chr)
x <- x[x$ground.target == chr, ]
d <- pivot_longer(x, chr13:chr22, "chromosome")

p <- ggplot(d, aes(x=start, y=value, color=chromosome)) +
  geom_step() +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0, 0, 0)) +
  scale_y_continuous(limits = c(0, max(10, max(d$value))), breaks=pretty_breaks()) +
  facet_wrap(~chromosome, scales = "free", ncol=1, labeller = labeller(variable = labels)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),

    text = element_text(size = 20),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),

    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.position = "top",

    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    plot.margin = unit(c(0,1,0,5), "cm")
  ) + labs(
  x = paste(chr, 'position'),
  y = paste('# contigs\n')
) +
  guides(colour = guide_legend(override.aes = list(size=10)))
#p
p_all_together <- ggplot(d, aes(x=start, y=value, color=chromosome)) +
  geom_step(alpha=0.5) +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0, 0, 0)) +
  scale_y_continuous(limits = c(0, max(10, max(d$value))), breaks=pretty_breaks()) +
  #facet_wrap(~chromosome, scales = "free", ncol=1, labeller = labeller(variable = labels)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    
    text = element_text(size = 20),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.position = "top",
    
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    plot.margin = unit(c(0,1,0,5), "cm")
  ) + labs(
    x = paste(chr, 'position'),
    y = paste('# contigs\n')
  ) +
  guides(colour = guide_legend(override.aes = list(size=10)))
#p_all_together

library(png)
library(grid)
img <- readPNG(path_annotation)

ggplotted_img <- ggplot() +
  annotation_custom(
    rasterGrob(img, width = 1, height = 1),
    xmin = - Inf, xmax = Inf,
    ymin = - Inf, ymax = Inf
  ) + theme(

  plot.margin = unit(c(0,1,0,0), "cm")
)

library(ggpubr)
p_with_annotation <- ggpubr::ggarrange(
  ggplotted_img, p,
  labels=c('', ''),
  heights = c(1.8, 3),
  legend = "right", # legend position,
  common.legend = T,
  nrow = 2
)
ggsave(
  plot = p_with_annotation,
  paste0(prefix_output, '.separated.pdf'),
  width = width, height = 7 * 4,
  units = "cm",
  dpi = 100, bg = "transparent",
  limitsize = FALSE
)


p_all_together_with_annotation <- ggpubr::ggarrange(
  ggplotted_img, p_all_together,
  labels=c('', ''),
  heights = c(1.8, 1),
  legend = "right", # legend position,
  common.legend = T,
  nrow = 2
)
ggsave(
  plot = p_all_together_with_annotation,
  paste0(prefix_output, '.together.pdf'),
  width = width, height = 7 * 2,
  units = "cm",
  dpi = 100, bg = "transparent",
  limitsize = FALSE
)
