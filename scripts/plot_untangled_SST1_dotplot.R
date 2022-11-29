args <- commandArgs()
path_untangle_tsv <- args[6]
label_sst1 <- args[7]
x_min_sst1 <- as.numeric(args[8])
x_max_sst1 <- as.numeric(args[9])
size_range <- as.numeric(args[10])
nth_best <- as.numeric(args[11])
path_query_to_consider <- args[12]
dir_output <- args[13]

library(ggplot2)
library(ggforce)
library(tidyverse)

query_to_consider <- read.delim(path_query_to_consider, header = F)

x_min <- x_min_sst1-1000000
x_max <- x_max_sst1+1000000

x <- read.delim(path_untangle_tsv) %>%
  rename(query.name=X.query.name) %>%
  filter(query.name %in% query_to_consider$V1 & nth.best <= nth_best & ref.start >= x_min & ref.end <= x_max)

#x <- x[x$self.coverage <= 1,]

# From https://doi.org/10.1093/bioinformatics/btac244
#x$estimated_identity <- exp((1.0 + log(2.0 * x$jaccard/(1.0+x$jaccard)))-1.0)

#x <- x[x$estimated_identity >= estimated_identity_threshold,]

options(scipen = 9)

ref.name <- unique(x$ref.name)

dir.create(dir_output, recursive=T)

for(q in unique(x$query.name)){
  for(nth_b in seq(1, nth_best)) {
    print(paste(q, nth_b))
    
    yy <- x %>% filter(query.name %in% q & nth.best <= nth_b)
    
    # Swap columns for inverted segments
    yy$ref.start_ = ifelse(yy$inv == '-', yy$ref.end, yy$ref.start)
    yy$ref.end_ = ifelse(yy$inv == '-', yy$ref.start, yy$ref.end)
    
    p <- ggplot(
      yy,
      aes(x = query.start, xend = query.end, y = ref.start_, yend = ref.end_)) +
      geom_segment(size = 0.3) +
      #facet_wrap(ref.name ~ query.name) +
      coord_fixed() +
      theme_bw() +
      theme(
        text = element_text(size = 12.6),
        axis.text.x = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 12),
        
        plot.title = element_text(hjust = 0.5)
        
        #panel.grid.minor = element_line(size = 0.125),
        #panel.grid.major = element_line(size = 0.25)
      ) +
      xlab(q) +
      ylab(ref.name) +
      ggtitle(paste0(q, ' vs ', ref.name, '; ', nth_b, ' best hit(s)')) +
      scale_y_reverse()
    p <- p + annotate("rect",
                      xmin = -Inf, xmax = Inf,
                      ymin = x_min_sst1, ymax = x_max_sst1,
                      fill = "#EE2222", alpha = .2, color = "#444444", size = 0.1)
    ann_text <- data.frame(
      `query.start` = min(yy$query.start)+(max(yy$query.start)-min(yy$query.start))/2,
      `query.end` = min(yy$query.start)+(max(yy$query.start)-min(yy$query.start))/2,
      `ref.start_` = x_min_sst1-(x_max_sst1-x_min_sst1),
      `ref.end_` = x_min_sst1-(x_max_sst1-x_min_sst1),
      `query.name` = factor(q, levels = q)
    )
    p <- p + geom_text(data = ann_text, label = label_sst1, size = 4)
    ggsave(plot = p, file.path(dir_output, paste0(q, '.n', nth_b, '.pdf')), width = 20, height = 20, units = "cm", dpi = 300, bg = "transparent", limitsize = FALSE)
  }
}
