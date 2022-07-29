args <- commandArgs()
path_entropy_tsv <- args[6]
x_min <- as.numeric(args[7])
x_max <- as.numeric(args[8])
width <- as.numeric(args[9])
num_chr <- as.numeric(args[10])
path_annotation <- args[11]
path_output <- args[12]

library(ggplot2)
library(ggforce)
library(scales) # for pretty_breaks()


x <- read.delim(path_entropy_tsv)

x[x$shannon_div_index == -1, ]$shannon_div_index <- NA

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

library(dplyr)
# Compute average SDI by window by ground.target
#y <- x %>% 
#  dplyr::group_by(ground.target, start.pos, end.pos) %>% 
#  dplyr::summarize(average_sdi = mean(shannon_div_index, na.rm = TRUE))

# Apply filters
chr <- paste0('chm13#chr', num_chr)
xx <- x[x$ground.target == chr,]


p <- ggplot(xx, aes(x=start, y=shannon_div_index, color=ground.target)) +
  geom_step(aes(alpha=num.queries)) +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0, 0, 0)) +
  scale_y_continuous(limits = c(0, max(2, max(xx$shannon_div_index))), breaks=pretty_breaks(n=6)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    
    text = element_text(size = 26),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.position = "top",
    
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    plot.margin = unit(c(0,1.03,0,7.1), "cm")
  ) + labs(
    x = paste('Position'),
    y = paste('Order entropy\n'),
    color = 'Target',
    alpha = '# contigs'
  ) +
  scale_color_manual(values=colors) +
  guides(colour = guide_legend(override.aes = list(size=10)))
#p


library(png)
library(grid)
img <- readPNG(path_annotation)

ggplotted_img <- ggplot() +
  annotation_custom(
    rasterGrob(img, width = 1, height = 1),
    xmin = - Inf, xmax = Inf,
    ymin = - Inf, ymax = Inf
  ) + theme(
    plot.margin = unit(c(0,0.95,0.5,0.3), "cm")
  )

library(ggpubr)
p_with_annotation <- ggpubr::ggarrange(
  ggplotted_img, p,
  labels=c('', ''),
  heights = c(8, max(7, length(unique(xx$contig)))),
  legend = "right", # legend position,
  common.legend = T,
  nrow = 2
)
ggsave(
  plot = p_with_annotation,
  path_output,
  width = width, height = (6+max(7, length(unique(xx$ground.target)))),
  units = "cm",
  dpi = 100, bg = "white",
  limitsize = FALSE
)

