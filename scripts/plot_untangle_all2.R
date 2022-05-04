args <- commandArgs()
path_untangle_grounded_all_tsv <- args[6]
title <- args[7]
x_min <- as.numeric(args[8])
x_max <- as.numeric(args[9])
width <- as.numeric(args[10])
height <- as.numeric(args[11])
nth.best <- as.numeric(args[12])
ref.nth.best <- as.numeric(args[13])

library(ggplot2)
library(ggforce)

x <- read.delim(path_untangle_grounded_all_tsv)

# To have it as numeric column
#x$self.coverage[x$self.coverage == '.'] <- 1
#x$self.coverage <- as.numeric(x$self.coverage)

#x <- x[x$self.coverage <= 1,]
x <- x[x$nth.best  <= nth.best & x$ref.nth.best <= ref.nth.best,]

colors <- c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF4")
if (length(unique(x$target)) > 5) {
  # Mobin's annotations
  colors <- c(colors, "#000000", "#000000", "#000000")
}

p <- ggplot(
    x,
    aes(
      x = ref.begin + (ref.end - ref.begin) / 2, width = ref.end - ref.begin - 200, y = query, fill = target#, alpha=self.coverage
    )
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
  xlim(x_min, x_max) +
  scale_fill_manual(values=colors)

ggsave(plot = p, paste0(path_untangle_grounded_all_tsv, '.pdf'), width = width, height = height, units = "cm", dpi = 100, bg = "transparent", limitsize = FALSE)
