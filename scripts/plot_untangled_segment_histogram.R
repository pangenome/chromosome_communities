args <- commandArgs()
path_untangle_grounded_tsv <- args[6]
x_min <- as.numeric(args[7])
x_max <- as.numeric(args[8])
width <- as.numeric(args[9])
height <- as.numeric(args[10])
max_len <- as.numeric(args[11])
nth.best <- as.numeric(args[12])
ref.nth.best <- as.numeric(args[13])
estimated_identity_threshold <- as.numeric(args[14])
path_query_to_consider <- args[15]
path_output <- args[16]


library(ggplot2)
library(ggforce)
library(tidyverse)

query_to_consider <- read.delim(path_query_to_consider, header = F)

x <- read.delim(path_untangle_grounded_tsv) %>%
  filter(query %in% query_to_consider$V1)

# To avoid errors
if (sum(x$jaccard > 1) > 0) {
  x[x$jaccard > 1,]$jaccard <- 1
}

# From https://doi.org/10.1093/bioinformatics/btac244
x$estimated_identity <- exp((1.0 + log(2.0 * x$jaccard/(1.0+x$jaccard)))-1.0)

x <- x[x$estimated_identity >= estimated_identity_threshold,]

colors <- c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF4")
if (length(unique(x$target)) > 5) {
  # Mobin's annotations
  colors <- c(colors, "#000000", "#000000", "#000000")
}

# Apply filters
xx <- x[x$nth.best <= nth.best & x$ref.nth.best <= ref.nth.best,]

# Do not consider dedicated annotation bars
# Do not consider other acros references
xx <- xx %>%
  filter(!grepl('rDNA', query) & !grepl('centromere', query)) %>%
  filter(!grepl('chr', query))

# To group by query
xx$query.hacked <- paste(xx$query, xx$nth.best, sep = "-")

xx <- xx %>%
  arrange(query.hacked)

# Before as.character, then as.numeric, else the numbers will correspond each factor level (1, 2) rather than the levels themselves!
xx$ref.end <- as.numeric(as.character(xx$ref.end))
xx$ref.begin <- as.numeric(as.character(xx$ref.begin))
xx$ref_len <- xx$ref.end - xx$ref.begin

xx <- xx[xx$ref.end >= x_min && xx$ref.end <= x_max,]

p <- ggplot(xx, aes(x=ref_len)) + 
  geom_histogram(aes(fill=grounded.target)) +
  #geom_density(alpha=.2, fill="#FF6666") +
  facet_grid(~grounded.target, scales = "free_y", space = "free") +
  theme_bw() +
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
  scale_x_continuous(limits = c(x_min, max_len), expand = c(0, 0)) +
  scale_fill_manual(name = "Target", values = colors) +
  labs(x = "Length")

ggsave(plot = p, path_output, width = width, height = height, units = "cm", dpi = 100, bg = "transparent", limitsize = FALSE)
