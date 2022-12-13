args <- commandArgs()
path_untangle_tsv <- args[6]
label_sst1 <- args[7]
x_min_sst1 <- as.numeric(args[8])
x_max_sst1 <- as.numeric(args[9])
size_range <- as.numeric(args[10])
nth_best <- as.numeric(args[11])
path_contig_to_len <- args[12]
#path_query_to_flip <- args[13]
path_query_to_consider <- args[13]
dir_output <- args[14]

library(ggplot2)
library(ggforce)
library(tidyverse)

query_to_consider <- read.delim(path_query_to_consider, header = F)

#contig_to_flip <- read.delim(path_contig_to_flip, header = F) %>%
#  rename(contig.name=V1)

contig2len <- read.delim(path_contig_to_len, header = F) %>%
  rename(contig.name=V1, len=V2)

x <- read.delim(path_untangle_tsv) %>%
  rename(query.name=X.query.name) %>%
  filter(query.name %in% query_to_consider$V1)

x_min <- x_min_sst1-size_range
x_max <- x_max_sst1+size_range

# Flip coordinates
for(q in unique(x$query.name)){
  #if (q %in% contig_to_flip$contig.name){
    print(q)
    
    # Take contig length
    len <- contig2len %>% filter(contig.name == q) %>% pull(len)
    
    # Flip coordinates
    tmp <- x[x$query.name == q,]$query.start
    x[x$query.name == q,]$query.start <- len - x[x$query.name == q,]$query.end
    x[x$query.name == q,]$query.end <- len - tmp
  #}
}
for(r in unique(x$ref.name)){
  print(r)
  
  #if (r %in% contig_to_flip$contig.name){
    # Take contig length
    len <- contig2len %>% filter(contig.name == r) %>% pull(len)
    
    # Flip coordinates
    tmp <- x[x$ref.name == r,]$ref.start
    x[x$ref.name == r,]$ref.start <- len - x[x$ref.name == r,]$ref.end
    x[x$ref.name == r,]$ref.end <- len - tmp
  #}
}

xx <- x %>% filter(nth.best <= nth_best & ref.start >= x_min & ref.end <= x_max)

#x <- x[x$self.coverage <= 1,]

# From https://doi.org/10.1093/bioinformatics/btac244
xx$estimated_identity <- exp((1.0 + log(2.0 * xx$score/(1.0+xx$score)))-1.0)

#x <- x[x$estimated_identity >= estimated_identity_threshold,]

options(scipen = 9)

ref.name <- unique(x$ref.name)

dir.create(dir_output, recursive=T)

for(q in unique(x$query.name)){
  for(nth_b in c(1, nth_best)) {
    print(paste(q, nth_b))

    yy <- xx %>% filter(query.name %in% q & nth.best <= nth_b)
    
    # Swap columns for inverted segments
    yy$ref.start_ = ifelse(yy$inv == '-', yy$ref.end, yy$ref.start)
    yy$ref.end_ = ifelse(yy$inv == '-', yy$ref.start, yy$ref.end)
    
    p <- ggplot(
      yy,
      aes(x = query.start, xend = query.end, y = ref.start_, yend = ref.end_, color = estimated_identity)) +
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
      labs(color='Estimated identity') +
      scale_y_reverse() + scale_x_continuous(position = "top")
    p <- p + annotate("rect",
                      xmin = -Inf, xmax = Inf,
                      ymin = x_min_sst1, ymax = x_max_sst1,
                      fill = "#EE2222", alpha = .2, color = "#444444", size = 0.1)
    ann_text <- data.frame(
      `query.start` = min(yy$query.start)+(max(yy$query.start)-min(yy$query.start))/2,
      `query.end` = min(yy$query.start)+(max(yy$query.start)-min(yy$query.start))/2,
      `ref.start_` = x_min_sst1+(x_max_sst1-x_min_sst1)/2,
      `ref.end_` = x_min_sst1+(x_max_sst1-x_min_sst1)/2,
      `query.name` = factor(q, levels = q),
      `estimated_identity` = 1.0
    )
    p <- p + geom_text(data = ann_text, label = label_sst1, size = 4)

    ggsave(plot = p, file.path(dir_output, paste0(q, '.vs.', ref.name, '.n', nth_b, '.pdf')), width = 20, height = 20, units = "cm", dpi = 300, bg = "transparent", limitsize = FALSE)
  }
}
