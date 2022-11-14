args <- commandArgs()
path_untangle_grounded_all_tsv <- args[6]
title <- args[7]
x_min <- as.numeric(args[8])
x_max <- as.numeric(args[9])
width <- as.numeric(args[10])
height <- as.numeric(args[11])

library(ggplot2)
library(ggforce)
library(tidyverse)

panel_spacing <- 0.1
nth.best <- 1
ref.nth.best <- 1

x <- read.delim(path_untangle_grounded_all_tsv)

xx <- x[x$nth.best <= nth.best & x$ref.nth.best <= ref.nth.best,]
# To have it as numeric column
#xx$self.coverage[xx$self.coverage == '.'] <- 1
#xx$self.coverage <- as.numeric(xx$self.coverage)

#xx <- xx[xx$self.coverage <= 1,]

# Do not consider dedicated annotation bars
# Do not consider other acros references
xx <- xx %>%
  filter(!grepl('rDNA', query) & !grepl('centromere', query)) %>%
  filter(!grepl('chr', query)) %>%
  filter(!grepl('HG002', query))

# To group by query
xx$query.hacked <- paste(xx$query, xx$nth.best, sep = "-")
xx$query2 <- gsub('_','',xx$query)

xx <- xx %>%
  arrange(query.hacked)

colors <- c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF4")
if (length(unique(xx$target)) > 5) {
  # Mobin's annotations
  colors <- c(colors, "#000000", "#000000", "#000000", "#000000")
}

p <- ggplot(
    xx,
    aes(
      x = ref.begin + (ref.end - ref.begin) / 2, width = ref.end - ref.begin,
      y = ordered(query, levels = rev(unique(query))),
      fill = target#, alpha=self.coverage
    )
  ) +
  geom_tile() +
  ggtitle(title) +
  facet_col(query2~grounded.target, scales = "free_y", space = "free", labeller = labeller(variable = labels)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    
    text = element_text(size = 22),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22),
    legend.position = "top",
    
    panel.spacing = unit(panel_spacing, "lines"),
    #panel.border = element_rect(color = "grey", fill = NA, size = 1), #element_blank(),
    
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.y = element_blank(),
    
    plot.margin = unit(c(0,1.03,0,2.20), "cm"),
  ) +
  xlim(x_min, x_max) +
  scale_fill_manual(values=colors)

ggsave(plot = p, paste0(path_untangle_grounded_all_tsv, '.pdf'), width = width, height = height, units = "cm", dpi = 100, bg = "transparent", limitsize = FALSE)
