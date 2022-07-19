args <- commandArgs()
path_untangle_bed <- args[6]
path_output <- args[7]

library(ggplot2)
#library(knitr)
x <- read.delim(path_untangle_bed, sep = '\t', header = T, comment.char = "$")
colnames(x)[1] <- 'query.name'

#x$query.name <- gsub(":.*", "", x$query.name)
#x$query.name <- gsub("#J.*", "", x$query.name)

x <- x[x$ref.start < 25000000 & x$ref.end <=25000000,]

# To avoid errors
if (sum(x$jaccard > 1) > 0) {
  x[x$jaccard > 1,]$jaccard <- 1
}

# From https://doi.org/10.1093/bioinformatics/btac244
x$estimated_identity <- exp((1.0 + log(2.0 * x$jaccard/(1.0+x$jaccard)))-1.0)

x <- x[x$estimated_identity >= estimated_identity_threshold,]


#x_subset <- subset(x, query.name %in% c(
#       "HG002#MAT#chr13.prox"
#   )
#)

p <- ggplot(
  x, aes(x = query.start, xend = query.end, y = ref.start, yend = ref.end, color=ref.name)) +
  geom_segment(size = 1.2) +
  facet_grid(query.name ~ ref.name) +
  coord_fixed() +
  theme_bw() +
  theme(
    text = element_text(size = 12.6),
    axis.text.x = element_text(size = 12, angle = 90),
    axis.text.y = element_text(size = 12),
    
    #panel.grid.minor = element_line(size = 0.125),
    #panel.grid.major = element_line(size = 0.25)
  ) +
  xlab("Query start") +
  ylab("Reference start")  + theme(aspect.ratio=1)

ggsave(
    plot = p,
    file.path(path_output, 'DotPlots.0_25Mbps.pdf'),
    width = 70, height = 1400, units = "cm", dpi = 300, bg = "transparent", limitsize = F)
