args <- commandArgs()
path_untangle_grounded_all_tsv <- args[6]
x_min <- as.numeric(args[7])
x_max <- as.numeric(args[8])
width <- as.numeric(args[9])
nth.best <- as.numeric(args[10])
ref.nth.best <- as.numeric(args[11])
num_chr <- as.numeric(args[12])
path_annotation <- args[13]
path_output <- args[14]


library(ggplot2)
library(ggforce)
library(tidyverse)

x <- read.delim(path_untangle_grounded_all_tsv)

# To have it as numeric column
#x$self.coverage[x$self.coverage == '.'] <- 1
#x$self.coverage <- as.numeric(x$self.coverage)

#x <- x[x$self.coverage <= 1,]

colors <- c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF4")
if (length(unique(x$target)) > 5) {
  # Mobin's annotations
  colors <- c(colors, "#000000", "#000000", "#000000")
}

# Apply filters
chr <- paste0('chm13#chr', num_chr)
xx <- x[x$grounded.target == chr & x$nth.best <= nth.best & x$ref.nth.best <= ref.nth.best,]

# Do not display HG002 contigs
# Do not display dedicated annotation bars
# Do not display other acros references
xx <- xx %>%
  filter(!grepl('HG002#1', query) & !grepl('HG002#2', query)) %>%
  filter(!grepl('rDNA', query) & !grepl('centromere', query)) %>%
  filter(!grepl('chr', query) | grepl(paste0('chr', num_chr), query))

# To group by query
xx$query.hacked <- paste(xx$query, xx$nth.best, sep ='-')

xx <- xx %>%
  arrange(query.hacked)

p <- ggplot(
  xx,
  aes(
    x = ref.begin + (ref.end - ref.begin) / 2, width = ref.end - ref.begin - 200,
    y = ordered(query.hacked, levels=rev(unique(query.hacked))),
    fill = target,
    alpha=jaccard
  )
) +
  geom_tile() +
  #ggtitle(paste(chr, title)) +
  facet_grid(query~grounded.target, scales = "free_y", space = "free", labeller = labeller(variable = labels)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),

    text = element_text(size = 32),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),

    legend.title = element_text(size = 32),
    legend.text = element_text(size = 32),
    legend.position = "top",

    panel.spacing = unit(0.8, "lines"),
    #panel.border = element_rect(color = "grey", fill = NA, size = 1), #element_blank(),

    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.y=element_blank()
  ) +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
  scale_fill_manual(values=colors) +
  labs(x ="Position")
#ggsave(plot = p, paste0(path_untangle_grounded_all_tsv, '.pdf'), width = width, height = height, units = "cm", dpi = 100, bg = "transparent", limitsize = FALSE)


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
  heights = c(1, length(unique(xx$query))/3),
  legend = "right", # legend position,
  common.legend = T,
  nrow = 2
)

ggsave(
  plot = p_with_annotation,
  path_output,
  width = width, height = length(unique(xx$query)) * 4,
  units = "cm",
  dpi = 100, bg = "transparent",
  limitsize = FALSE
)
