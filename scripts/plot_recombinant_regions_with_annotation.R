args <- commandArgs()
path_recombinant_regions_table_tsv <- args[6]
x_min <- as.numeric(args[7])
x_max <- as.numeric(args[8])
width <- as.numeric(args[9])
height <- as.numeric(args[10])
num_chr <- as.numeric(args[11])
path_annotation <- args[12]
path_output <- args[13]

library(ggplot2)

xx <- read.delim(path_recombinant_regions_table_tsv, header = F)
colnames(xx) <- c('self.coverage.threshold', 'estimated.identity.threshold',  'grounded.target', 'ref.begin', 'ref.end', 'num.contigs')

# Remove ranges not supported by any contigs
xx <- xx[xx$num.contigs > 0,]

xx$label <- paste0('sc', xx$self.coverage.threshold, '.eid', xx$estimated.identity.threshold)

# Total size of each region
#library(tidyverse)
#xx$len <- xx$ref.end-xx$ref.begin
#xx %>% 
#  group_by(label) %>% 
#  summarise(Total = sum(len)) %>% View()

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


p <- ggplot(
  xx,
  aes(
    x = ref.begin + (ref.end - ref.begin) / 2, width = ref.end - ref.begin,
    y = label,
    fill = grounded.target,
    alpha = num.contigs
    #y = ordered(query, levels = rev(unique(query))),
  )
) +
  geom_tile() +
  #ggtitle(paste(chr, title)) +
  facet_grid(label~grounded.target, scales = "free_y", space = "free", labeller = labeller(variable = labels)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    
    text = element_text(size = 32),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    
    legend.title = element_text(size = 32),
    legend.text = element_text(size = 32),
    legend.position = "top",
    
    panel.spacing = unit(0.0, "lines"),
    #panel.border = element_rect(color = "grey", fill = NA, size = 1), #element_blank(),
    
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.y = element_blank(),
    
    plot.margin = unit(c(0,1.03,0,6.7), "cm"),
  ) +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
  scale_fill_manual(values = colors) +
  #scale_alpha_continuous(limits = c(1,31), breaks = c(0,25,50,100)) +
  scale_alpha_continuous(breaks = round(seq(1, max(xx$num.contigs), by = 5), 3)) +
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
  heights = c(1, 1.7),
  legend = "right", # legend position,
  common.legend = T,
  nrow = 2
)

ggsave(
  plot = p_with_annotation,
  path_output,
  width = width, height = height,
  units = "cm",
  dpi = 100, bg = "white",
  limitsize = FALSE
)

