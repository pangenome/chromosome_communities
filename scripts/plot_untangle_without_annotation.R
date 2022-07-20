args <- commandArgs()
path_untangle_grounded_tsv <- args[6]
x_min <- as.numeric(args[7])
x_max <- as.numeric(args[8])
width <- as.numeric(args[9])
height_row <- as.numeric(args[10])
panel_spacing <- as.numeric(args[11])
nth.best <- as.numeric(args[12])
ref.nth.best <- as.numeric(args[13])
num_chr <- args[14]
estimated_identity_threshold <- as.numeric(args[15])
path_query_to_consider <- args[16]
path_output <- args[17]


library(ggplot2)
library(ggforce)
library(tidyverse)

query_to_consider <- read.delim(path_query_to_consider, header = F)

x <- read.delim(path_untangle_grounded_tsv) %>%
  filter(query %in% query_to_consider$V1)

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

if (length(unique(x$target)) == 2) {
  # Sex chromosomes
  colors <- c("#E76BF3", "#00BFC4")
} else {
  # Acrocentric chromosomes
  colors <- c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF4")
  if (length(unique(x$target)) > 5) {
    # Mobin's annotations
    colors <- c(colors, "#000000", "#000000", "#000000")
  }
}


# Apply filters
chr <- paste0('chm13#chr', num_chr)
xx <- x[x$grounded.target == chr & x$nth.best <= nth.best & x$ref.nth.best <= ref.nth.best,]

# Do not consider dedicated annotation bars
# Do not consider other acros references
xx <- xx %>%
  filter(!grepl('rDNA', query) & !grepl('centromere', query)) %>%
  filter(!grepl('chr', query) | grepl(paste0('chr', num_chr), query))

# To group by query
xx$query.hacked <- paste(xx$query, xx$nth.best, sep = "-")

xx <- xx %>%
  arrange(query.hacked)

# To avoid errors
if (sum(xx$jaccard > 1) > 0) {
  xx[xx$jaccard > 1,]$jaccard <- 1
}

# Filter before to avoid plotting empty rows
xx$x <- xx$ref.begin + (xx$ref.end - xx$ref.begin) / 2
xx <- xx[xx$x >= x_min & xx$x <= x_max,]

p <- ggplot(
  xx,
  aes(
    x = x, width = ref.end - ref.begin ,
    y = ordered(query.hacked, levels = rev(unique(query.hacked))),
    fill = target,
    alpha = estimated_identity
  )
) +
  geom_tile() +
  #ggtitle(paste(chr, title)) +
  facet_grid(query~grounded.target, scales = "free_y", space = "free", labeller = labeller(variable = labels)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),

    text = element_text(size = 32),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),

    legend.title = element_text(size = 32),
    legend.text = element_text(size = 32),
    legend.position = "top",

    panel.spacing = unit(panel_spacing, "lines"),
    #panel.border = element_rect(color = "grey", fill = NA, size = 1), #element_blank(),

    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
  scale_fill_manual(values = colors) +
  labs(x = "Position")
ggsave(plot = p, path_output, width = width, height = height_row*length(unique(xx$query))*nth.best, units = "cm", dpi = 100, bg = "transparent", limitsize = FALSE)
