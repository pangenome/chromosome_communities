args <- commandArgs()
path_entropy_tsv <- args[6]
x_min <- as.numeric(args[7])
x_max <- as.numeric(args[8])
width <- as.numeric(args[9])
height_bar <- as.numeric(args[10])
num_chr <- as.numeric(args[11])
path_annotation <- args[12]
path_output <- args[13]

library(ggplot2)
library(ggforce)

x <- read.delim(path_entropy_tsv)

x <- x[x$shannon_div_index != -1, ]

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


# Apply filters
chr <- paste0('chm13#chr', num_chr)
xx <- x[x$grounded.target == chr,]

p <- ggplot(xx,
    aes(
        x = start.pos + (end.pos - start.pos) / 2,
        width = end.pos - start.pos - 200,
        y = contig, alpha=shannon_div_index,
        fill = grounded.target
    )
   ) +
  geom_tile() +
  #ggtitle(title) +
  facet_col(~grounded.target, scales = "free_y", space = "free", labeller = labeller(variable = labels)) +
  theme(
    plot.title = element_text(hjust = 0.5),

    text = element_text(size = 32),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),

    legend.title = element_text(size = 32),
    legend.text = element_text(size = 32),
    legend.position = "top",

    #axis.title.y=element_blank()
  ) +
  labs(alpha="SDI", fill="Target") +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
  scale_fill_manual(values=colors) +
  labs(x ="Position")

#ggsave(plot = p, paste0(path_entropy_tsv, '.entropy.pdf'), width = 120, height = height, units = "cm", dpi = 100, bg = "transparent", limitsize = FALSE)
p


library(png)
library(grid)
img <- readPNG(path_annotation)

ggplotted_img <- ggplot() +
  annotation_custom(
    rasterGrob(img, width = 1, height = 1),
    xmin = - Inf, xmax = Inf,
    ymin = - Inf, ymax = Inf
  )

library(ggpubr)
p_with_annotation <- ggpubr::ggarrange(
  ggplotted_img, p,
  labels=c('', ''),
  heights = c(height_bar*8, height_bar*max(7, length(unique(xx$contig)))),
  legend = "right", # legend position,
  common.legend = T,
  nrow = 2
)

ggsave(
  plot = p_with_annotation,
  path_output,
  width = width, height = (8+max(7, length(unique(xx$contig)))) * height_bar,
  units = "cm",
  dpi = 100, bg = "white",
  limitsize = FALSE
)
