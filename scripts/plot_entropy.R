args <- commandArgs()
path_sdi_tsv <- args[6]
title <- args[7]

library(ggplot2)
library(ggforce)

x <- read.delim(path_sdi_tsv)

for(ground in unique(x$ground.target)) {
  x_sub <- subset(x, ground.target == ground)

  print(ground)

  p <- ggplot(
      x_sub,
      aes(x = start + (end - start) / 2, width = end - start - 200, y = shannon_div_index)
    ) +
    facet_col(~query, scales = "free_y", space = "free", labeller = labeller(variable = labels)) +
    geom_point() +
    ylim(0,max(x_sub$shannon_div_index) + 0.1) +
    ggtitle(title)

    num_queries <- length(unique(x_sub$query))
    ggsave(plot = p, paste0(path_sdi_tsv, '.', ground, '.pdf'), width = 80, height = num_queries * 8, units = "cm", dpi = 300, bg = "transparent", limitsize = F)
}

#query	ground.target	start	end	shannon_div_index
#p <- ggplot(x, aes(x = start + (end - start) / 2, width = end - start - 200, y = shannon_div_index)) +
#  facet_grid(query~ground.target) +
#  geom_point() +
#  ggtitle(title)

#num_queries <- length(unique(x$query))
#ggsave(plot = p, paste0(path_sdi_tsv, '.pdf'), width = 80, height = num_queries * 6, units = "cm", dpi = 300, bg = "transparent", limitsize = F)
