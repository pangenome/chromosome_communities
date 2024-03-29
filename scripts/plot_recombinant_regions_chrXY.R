args <- commandArgs()
path_chrXY_plus_recombinant <- args[6]
width <- as.numeric(args[7])
height <- as.numeric(args[8])
num_chr <- args[9]
path_output <- args[10]


library(ggplot2)
library(ggforce)
library(tidyverse)
library(scales) # for pretty_breaks()

x <- read.delim(path_chrXY_plus_recombinant)


if (num_chr == "X") {
  colors <- c("#E76BF3")
} else if (num_chr == "Y") {
  colors <- c("#00BFC4")
}else {
  colors <- c("#000000")
}

# Apply filters
chr <- paste0('chm13#chr', num_chr)
xx <- x[x$grounded.target == chr,]

# To group by query
xx$query.hacked <- xx$query

# Filter before to avoid plotting empty rows
xx$x <- xx$ref.begin + (xx$ref.end - xx$ref.begin) / 2
#xx <- xx[xx$x >= x_min & xx$x <= x_max,]

p <- ggplot(
  xx,
  aes(
    x = x, width = ref.end - ref.begin ,
    y = ordered(query.hacked, levels = rev(unique(query.hacked))),
    fill = target,
    alpha = self.coverage
  )
) +
  geom_tile() +
  #ggtitle(paste(chr, title)) +
 # facet_grid(~grounded.target, scales = "free_y", space = "free", labeller = labeller(variable = labels)) +
  theme_bw() + 
  scale_alpha(range=c(0.05,max(xx$self.coverage)/100.0), limits=c(0,max(xx$self.coverage)), breaks=seq(min(xx$self.coverage), 100, 20)) +
  theme(
    plot.title = element_text(hjust = 0.5),
    
    text = element_text(size = 32),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    
    legend.title = element_text(size = 32),
    legend.text = element_text(size = 32),
    legend.position = "top",
    
    #panel.spacing = unit(panel_spacing, "lines"),
    #panel.border = element_rect(color = "grey", fill = NA, size = 1), #element_blank(),
    
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  #scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
  scale_fill_manual(values = colors) +
  labs(x = "Position", alpha = '# contigs')
ggsave(plot = p, path_output, width = width, height = height, units = "cm", dpi = 100, bg = "transparent", limitsize = FALSE)

