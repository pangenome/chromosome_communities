args <- commandArgs()
path_entropy_tsv <- args[6]
title <- args[7]
x_max <- as.numeric(args[8])
height <- as.numeric(args[9])

library(ggplot2)
library(ggforce)

x <- read.delim(path_entropy_tsv)

x <- x[x$shannon_div_index != -1, ]

p <- ggplot(x, aes(x = start + (end - start) / 2, width = end - start, y = query, alpha=shannon_div_index, fill = ground.target)) +
  geom_tile() +
  ggtitle(title) +
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
  labs(alpha="Shannon Diversity Index", fill="Target") +
  xlim(0, x_max)

ggsave(plot = p, paste0(path_entropy_tsv, '.entropy.pdf'), width = 120, height = height, units = "cm", dpi = 100, bg = "transparent", limitsize = FALSE)
