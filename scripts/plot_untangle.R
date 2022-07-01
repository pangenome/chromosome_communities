args <- commandArgs()
path_untangle_grounded_tsv <- args[6]
title <- args[7]
x_max <- as.numeric(args[8])
height <- as.numeric(args[9])

library(ggplot2)

x <- read.delim(path_untangle_grounded_tsv)

colors <- c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF4")
if (length(unique(x$target)) > 5) {
  # Mobin's annotations
  colors <- c(colors, "#FF0000", "#222222", "#0000FF")
}

p <- ggplot(x, aes(x = ref.begin + (ref.end - ref.begin) / 2, width = ref.end - ref.begin, y = query, fill = target)) +
  geom_tile() +
  ggtitle(title) +
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
  xlim(0, x_max) +
  scale_fill_manual(values=colors)

ggsave(plot = p, paste0(path_untangle_grounded_tsv, '.pdf'), width = 120, height = height, units = "cm", dpi = 300, bg = "transparent", limitsize = FALSE)
