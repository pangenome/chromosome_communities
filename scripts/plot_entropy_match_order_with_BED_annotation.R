args <- commandArgs()
path_entropy_tsv <- args[6]
x_min <- as.numeric(args[7])
x_max <- as.numeric(args[8])
width <- as.numeric(args[9])
num_chr <- args[10]
path_annotation_bed <- args[11]
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
} else if (num_chr == 'X') {
  colors <- c("#E76BF3")
} else if (num_chr == 'Y') {
  colors <- c("#00BFC4")
}else {
  colors <- c("#000000")
}

library(dplyr)

# Apply filters
chr <- paste0('chm13#chr', num_chr)
xx <- x[x$ground.target == chr,]


p <- ggplot(xx, aes(x=start, y=shannon_div_index, color=ground.target)) +
  geom_step(aes(alpha=num.queries)) +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0, 0, 0)) +
  scale_y_continuous(limits = c(0, max(xx$shannon_div_index, na.rm = T) + 0.1), breaks=pretty_breaks(n=6)) +
  scale_alpha(limits=c(1,max(xx$num.queries))) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    
    text = element_text(size = 32),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.position = "top",
    
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    plot.margin = unit(c(0,0,0,0.5), "cm")
  ) + labs(
    x = paste('Position'),
    y = paste('Order entropy\n'),
    color = 'Target',
    alpha = '# contigs'
  ) +
  scale_color_manual(values=colors) +
  guides(colour = guide_legend(override.aes = list(size=10)))
#p


a <- read.delim(path_annotation_bed, header = F)
colnames(a) <- c('Target', 'ref.begin', 'ref.end')
a <- a[grepl(paste0('chr', num_chr), a$Target),]
a$ground.target <- chr

p_ann <- ggplot(
  a,
  aes(
    x = ref.begin + (ref.end - ref.begin) / 2,
    width = ref.end - ref.begin ,
    y = ordered(Target, levels = rev(unique(Target))),
    fill = ground.target
  )
) +
  geom_tile() +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    
    text = element_text(size = 32),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.position = "top",
    
    #panel.spacing = unit(panel_spacing, "lines"),
    #panel.border = element_rect(color = "grey", fill = NA, size = 1), #element_blank(),
    
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = unit(c(0,0,0.5,0.5), "cm")
  ) +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
  scale_fill_manual(values=colors) +
  labs(x = "Position")


library(ggpubr)
p_with_annotation <- ggpubr::ggarrange(
  p_ann, p,
  labels=c('', ''),
  heights = c(0.5, 1),
  legend = "right", # legend position,
  common.legend = T,
  nrow = 2
)
p_with_annotation
ggsave(
  plot = p_with_annotation,
  path_output,
  width = width, height = 10,
  units = "cm",
  dpi = 100, bg = "white",
  limitsize = FALSE
)
