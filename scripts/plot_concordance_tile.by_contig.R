args <- commandArgs()
path_concordance_tsv <- args[6]
title <- args[7]
x_max <- as.numeric(args[8])
height <- as.numeric(args[9])

library(ggplot2)
library(ggforce)

x <- read.delim(path_concordance_tsv)

x$concordance <- 2
x$concordance[x$contig.target == x$verkko.target] <- 1

x$concordance <- as.factor(x$concordance)

p <- ggplot(
  x,
  aes(
    x = start.pos + (end.pos - start.pos) / 2,
    width = end.pos - start.pos - 200,
    y = contig,
    alpha= concordance, fill = grounded.target)
  ) +
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
  labs(alpha="# different targets", fill="Target") +
  xlim(0, x_max)

ggsave(plot = p, paste0(path_concordance_tsv, '.concordance.pdf'), width = 120, height = height, units = "cm", dpi = 100, bg = "transparent", limitsize = FALSE)
