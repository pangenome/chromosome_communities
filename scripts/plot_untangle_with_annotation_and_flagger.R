args <- commandArgs()
path_untangle_grounded_tsv <- args[6]
x_min <- as.numeric(args[7])
x_max <- as.numeric(args[8])
width <- as.numeric(args[9])
height_bar <- as.numeric(args[10])
panel_spacing <- as.numeric(args[11])
nth.best <- as.numeric(args[12])
ref.nth.best <- as.numeric(args[13])
num_chr <- as.numeric(args[14])
estimated_identity_threshold <- as.numeric(args[15])
#path_query_to_consider <- args[16]
path_annotation <- args[16]
path_output <- args[17]


library(ggplot2)
library(ggforce)
library(tidyverse)

#query_to_consider <- read.delim(path_query_to_consider, header = F)

x <- read.delim(path_untangle_grounded_tsv)# %>%
  #filter(query %in% query_to_consider$V1)

# To have it as numeric column
#x$self.coverage[x$self.coverage == '.'] <- 1
#x$self.coverage <- as.numeric(x$self.coverage)

#x <- x[x$self.coverage <= 1,]

# To avoid errors
if (sum(x$jaccard > 1) > 0) {
  x[x$jaccard > 1,]$jaccard <- 1
}

# From https://doi.org/10.1093/bioinformatics/btac244
x$estimated_identity <- exp((1.0 + log(2.0 * x$jaccard/(1.0+x$jaccard)))-1.0)

x <- x[x$estimated_identity >= estimated_identity_threshold,]

colors <- c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF4")
if (length(unique(x$target)) > 5) {
  # Mobin's annotations
  x[x$target == 'Err',]$target = 'Unreliable'
  x[x$target == 'Unk',]$target = 'Unreliable'
  x[x$target == 'Col',]$target = 'Unreliable'
  x[x$target == 'Dup',]$target = 'Unreliable'
  x <- x %>% filter(target != 'Reliable')
  colors <- c(colors, '#262626')
}

# Apply filters
chr <- paste0('chm13#chr', num_chr)
xx <- x[x$grounded.target == chr & x$nth.best <= nth.best & x$ref.nth.best <= ref.nth.best,]

xx$grounded.target <- paste0('chr', num_chr)
  
# Do not consider dedicated annotation bars
# Do not consider other acros references
xx <- xx %>%
  filter(!grepl('rDNA', query) & !grepl('centromere', query)) %>%
  filter(!grepl('chr', query) & !grepl('HG002', query))

# To group by query
xx$query.hacked <- paste(xx$query, xx$nth.best, sep = "-")
xx$query2 <- gsub('_','',xx$query)

xx <- xx %>%
  arrange(query.hacked)

# To avoid errors
if (sum(xx$jaccard > 1) > 0) {
  xx[xx$jaccard > 1,]$jaccard <- 1
}

p <- ggplot(
  xx,
  aes(
    x = ref.begin + (ref.end - ref.begin) / 2, width = ref.end - ref.begin ,
    y = ordered(query, levels = rev(unique(query))),
    fill = target
  )
) +
  geom_tile() +
  #ggtitle(paste(chr, title)) +
  facet_grid(query2~grounded.target, scales = "free_y", space = "free", labeller = labeller(variable = labels)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),

    text = element_text(size = 22),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),

    legend.title = element_text(size = 32),
    legend.text = element_text(size = 32),
    legend.position = "top",

    panel.spacing = unit(panel_spacing, "lines"),
    #panel.border = element_rect(color = "grey", fill = NA, size = 1), #element_blank(),

    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.y = element_blank(),
    
    plot.margin = unit(c(0,1.03,0,2.20), "cm"),
  ) +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
  scale_fill_manual(values = colors) +
  labs(x = "Position", fill="Target")
p
#ggsave(plot = p, paste0(path_untangle_grounded_all_tsv, '.pdf'), width = width, height = height, units = "cm", dpi = 100, bg = "transparent", limitsize = FALSE)


library(png)
library(grid)
img <- readPNG(path_annotation)

ggplotted_img <- ggplot() +
  annotation_custom(
    rasterGrob(img, width = 1, height = 1),
    xmin = - Inf, xmax = Inf,
    ymin = - Inf, ymax = Inf
  ) + theme(
    plot.margin = unit(c(0,1,0.5,-0.3), "cm")
  )

library(ggpubr)
p_with_annotation <- ggpubr::ggarrange(
  ggplotted_img, p,
  labels=c('', ''),
  heights = c(height_bar*7.9, height_bar*length(unique(xx$query))*nth.best),
  legend = "right", # legend position,
  common.legend = T,
  nrow = 2
)

ggsave(
  plot = p_with_annotation,
  path_output,
  width = width, height = (12+length(unique(xx$query))*nth.best) * height_bar,
  units = "cm",
  dpi = 100, bg = "white",
  limitsize = FALSE
)
