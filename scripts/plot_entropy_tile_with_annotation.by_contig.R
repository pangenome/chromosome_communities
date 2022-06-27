args <- commandArgs()
path_entropy_tsv <- args[6]
x_min <- as.numeric(args[7])
x_max <- as.numeric(args[8])
width <- as.numeric(args[9])
height_bar <- as.numeric(args[10])
num_chr <- as.numeric(args[11])
path_annotation <- args[12]
prefix_output <- args[13]

library(ggplot2)
library(ggforce)

x <- read.delim(path_entropy_tsv)

x[x$shannon_div_index == -1, ]$shannon_div_index = NA

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
y <- x %>% 
  dplyr::group_by(ground.target, start.pos, end.pos) %>% 
  dplyr::summarize(average_sdi = mean(shannon_div_index, na.rm = TRUE))

# Apply filters
chr <- paste0('chm13#chr', num_chr)
xx <- x[x$ground.target == chr,]
yy <- y[y$ground.target == chr,]

p <- ggplot(xx[!is.na(xx$shannon_div_index),],
    aes(
        x = start.pos + (end.pos - start.pos) / 2,
        width = end.pos - start.pos - 200,
        y = contig, alpha=shannon_div_index,
        fill = ground.target
    )
   ) +
  geom_tile() +
  #ggtitle(title) +
  facet_col(~ground.target, scales = "free_y", space = "free", labeller = labeller(variable = labels)) +
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
#p


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
  paste0(prefix_output, '.by_contig.pdf'),
  width = width, height = (8+max(7, length(unique(xx$contig)))) * height_bar,
  units = "cm",
  dpi = 100, bg = "white",
  limitsize = FALSE
)


# Aggregate entropy by chromosome (averaging the SDI across contigs for each window)
if (FALSE) {
  p <- ggplot(yy[!is.na(yy$average_sdi),],
              aes(
                x = start.pos + (end.pos - start.pos) / 2,
                width = end.pos - start.pos - 200,
                y = ground.target, alpha=average_sdi,
                fill = ground.target
              )
  ) +
    geom_tile() +
    #ggtitle(title) +
    facet_col(~ground.target, scales = "free_y", space = "free", labeller = labeller(variable = labels)) +
    theme(
      plot.title = element_text(hjust = 0.5),
      
      text = element_text(size = 32),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      
      legend.title = element_text(size = 32),
      legend.text = element_text(size = 32),
      legend.position = "top",
      
      plot.margin = margin(0, 0.2, 0, 3.07, "cm")
      
      #axis.title.y=element_blank()
    ) +
    labs(alpha="Average SDI", fill="Target") +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    scale_fill_manual(values=colors) +
    labs(x ="Position", y = '')
  #p
}


library(scales) # for pretty_breaks()
p <- ggplot(yy, aes(x=start.pos + (end.pos - start.pos) / 2, y=average_sdi, color=ground.target)) +
  geom_step() +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0, 0, 0)) +
  scale_y_continuous(limits = c(0, max(2, max(y$average_sdi))), breaks=pretty_breaks()) +
  facet_wrap(~ground.target, scales = "free", ncol=1, labeller = labeller(variable = labels)) +
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
    plot.margin = unit(c(0,0.4,0,6.3), "cm")
  ) + labs(
    x = paste(chr, 'position'),
    y = paste('')
  ) +
  guides(colour = guide_legend(override.aes = list(size=10)))
#p

p_with_annotation <- ggpubr::ggarrange(
  ggplotted_img, p,
  labels=c('', ''),
  heights = c(height_bar*8, height_bar*max(7, length(unique(yy$ground.target))) / 1.7),
  legend = "right", # legend position,
  common.legend = T,
  nrow = 2
)
ggsave(
  plot = p_with_annotation,
  paste0(prefix_output, '.by_chromosome.pdf'),
  width = width, height = (6+max(7, length(unique(yy$ground.target)))) * height_bar,
  units = "cm",
  dpi = 100, bg = "white",
  limitsize = FALSE
)
