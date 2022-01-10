args <- commandArgs()
path_untangle_grounded_tsv <- args[6]
title <- args[7]

library(ggplot2)

x <- read.delim(path_untangle_grounded_tsv)

p <- ggplot(x, aes(x = ref.begin + (ref.end - ref.begin) / 2, width = ref.end - ref.begin - 200, y = query, fill = target)) +
  geom_tile() +
  ggtitle(title)

num_queries <- length(unique(x$query))
ggsave(plot = p, paste0(path_untangle_grounded_tsv, '.pdf'), width = 40, height = num_queries / 2, units = "cm", dpi = 300, bg = "transparent")
