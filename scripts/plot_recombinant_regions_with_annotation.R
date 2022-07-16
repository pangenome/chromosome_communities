library(ggplot2)

num_chr <- 13
path_annotation <- paste0('/home/guarracino/git/chromosome_communities/data/annotation/hgt_genome_euro_chr', num_chr, '_0_25Mbp.png')
path_output <- paste0('/home/guarracino/', num_chr, '.png')

xx <- read.delim('/home/guarracino/git/chromosome_communities/x.tsv', header = F)
colnames(xx) <- c('query', 'ref.begin', 'ref.end', 'grounded.target')

# Apply filters
chr <- paste0('chm13#chr', num_chr)
xx <- xx[xx$grounded.target == chr,]


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

x_min = 0
x_max = 25000000



p <- ggplot(
  xx,
  aes(
    x = ref.begin + (ref.end - ref.begin) / 2, width = ref.end - ref.begin ,
    y = ordered(query, levels = rev(unique(query))),
    fill = grounded.target,
    #alpha = jaccard
  )
) +
  geom_tile() +
  #ggtitle(paste(chr, title)) +
  #facet_grid(query~grounded.target, scales = "free_y", space = "free", labeller = labeller(variable = labels)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    
    text = element_text(size = 32),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    
    legend.title = element_text(size = 32),
    legend.text = element_text(size = 32),
    legend.position = "top",
    
    #panel.spacing = unit(panel_spacing, "lines"),
    #panel.border = element_rect(color = "grey", fill = NA, size = 1), #element_blank(),
    
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.y = element_blank(),
    
    plot.margin = unit(c(0,1.03,0,6.88), "cm"),
  ) +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
  scale_fill_manual(values = colors) +
  labs(x = "Position", fill="Target")



library(png)
library(grid)
img <- readPNG(path_annotation)

ggplotted_img <- ggplot() +
  annotation_custom(
    rasterGrob(img, width = 1, height = 1),
    xmin = - Inf, xmax = Inf,
    ymin = - Inf, ymax = Inf
  ) + theme(
    plot.margin = unit(c(0,1,0.5,0), "cm")
  )


library(ggpubr)
p_with_annotation <- ggpubr::ggarrange(
  ggplotted_img, p,
  labels=c('', ''),
  heights = c(3, 1),
  legend = "right", # legend position,
  common.legend = T,
  nrow = 2
)

width = 90
height = 10
ggsave(
  plot = p_with_annotation,
  path_output,
  width = width, height = height,
  units = "cm",
  dpi = 100, bg = "white",
  limitsize = FALSE
)
