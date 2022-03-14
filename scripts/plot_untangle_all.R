args <- commandArgs()
path_untangle_grounded_all_tsv <- args[6]
title <- args[7]
x_max <- args[8]
height <- args[9]

library(ggplot2)
library(ggforce)

x <- read.delim(path_untangle_grounded_all_tsv)

p <- ggplot(x, aes(x = ref.begin + (ref.end - ref.begin) / 2, width = ref.end - ref.begin - 200, y = query, fill = target)) +
  geom_tile() +
  ggtitle(title) +
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
  xlim(0, x_max)

ggsave(plot = p, paste0(path_untangle_grounded_all_tsv, '.pdf'), width = 120, height = height, units = "cm", dpi = 100, bg = "transparent", limitsize = FALSE)
