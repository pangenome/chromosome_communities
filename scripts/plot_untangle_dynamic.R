path_untangle_grounded_tsv <- '/home/guarracino/Downloads/Pangenomics/chromosome_communities/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.ALL.e50000.m1000.grounded.pq_touching.reliable.tsv.gz'
x_min <-  9000000
x_max <- 11000000
nth.best <- 3
ref.nth.best <- 1
num_chr <- 13
path_query_to_consider <- args[14]
dir_output <- args[15]


library(ggplot2)
library(ggforce)
library(tidyverse)
library(plotly)

query_to_consider <- read.delim(path_query_to_consider, header = F)

query_to_consider <- c(
  'chm13#chr13',
  'grch38#chr13',
  'HG002#MAT#chr13.prox',
  'HG002#PAT#chr13.prox',
  'HG01361#2#JAGYYW010000010.1'
  #'HG01978#1#JAGYVS010000056.1',
  #'HG02486#1#JAGYVM010000043.1',
  #'HG03540#2#JAGYVX010000153.1'
  )

x <- read.delim(path_untangle_grounded_tsv) %>%
  filter(query %in% query_to_consider$V1)

# To have it as numeric column
#x$self.coverage[x$self.coverage == '.'] <- 1
#x$self.coverage <- as.numeric(x$self.coverage)

#x <- x[x$self.coverage <= 1,]

colors <- c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF4")
if (length(unique(x$target)) > 5) {
  # Mobin's annotations
  colors <- c(colors, "#000000", "#000000", "#000000")
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

# Before as.character, then as.numeric, else the numbers will correspond each factor level (1, 2) rather than the levels themselves!
xx$query.end <- as.numeric(as.character(xx$query.end))
xx$query.begin <- as.numeric(as.character(xx$query.begin))
xx$query.len <- xx$query.end - xx$query.begin

xx$target.end <- as.numeric(as.character(xx$target.end))
xx$target.begin <- as.numeric(as.character(xx$target.begin))
xx$target.len <- xx$target.end - xx$target.begin

xx$ref.end <- as.numeric(as.character(xx$ref.end))
xx$ref.begin <- as.numeric(as.character(xx$ref.begin))
xx$ref.len <- xx$ref.end - xx$ref.begin


# All queries
xxx <- xx[xx$ref.begin >= x_min & xx$ref.end <= x_max,] %>% 
  mutate(Info = paste(
    "\n - query: ", query, ':', query.begin, "-", query.end, " (", query.len, ")",
    "\n - target: ", target, ':', target.begin, "-", target.end, " (", target.len, ")",
    "\n - jaccard: ", jaccard,
    "\n - ground: ", grounded.target, ":", ref.begin, "-", ref.end, " (", ref.len, ")",
    sep = "")
  ) %>%
  mutate(x = ref.begin + (ref.end - ref.begin) / 2) %>%
  mutate(width = ref.end - ref.begin) %>%
  mutate(y = ordered(query.hacked, levels = rev(unique(query.hacked))))

p <- ggplot(
  xxx,
  aes(
    x = x,
    y = y,
    width = width,
    fill = target,
    alpha = jaccard,
    label = Info
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
    
    panel.spacing = unit(0.05, "lines"),
    
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
  scale_fill_manual(values = colors) +
  labs(x = "Position")
p


l <- ggplotly(p, tooltip = c('label') )%>%
  layout(legend = list(orientation = "h", x = 0, y =-0.1))
l$height <- 400 + length(unique(xxx$query)) * nth.best * 8
#l

htmlwidgets::saveWidget(l, paste0(chr, ".html"))



# One query at a time
for (query in unique(xx$query)){
  print(query)
  
  xxx <- xx[xx$query == query & xx$ref.begin >= x_min & xx$ref.end <= x_max,] %>% 
    mutate(Info = paste(
      "\n - query: ", query, ':', query.begin, "-", query.end, " (", query.len, ")",
      "\n - target: ", target, ':', target.begin, "-", target.end, " (", target.len, ")",
      "\n - jaccard: ", jaccard,
      "\n - ground: ", grounded.target, ":", ref.begin, "-", ref.end, " (", ref.len, ")",
      sep = "")
    ) %>%
    mutate(x = ref.begin + (ref.end - ref.begin) / 2) %>%
    mutate(width = ref.end - ref.begin ) %>%
    mutate(y = ordered(query.hacked, levels = rev(unique(query.hacked))))
  
  p <- ggplot(
    xxx,
    aes(
      x = x,
      y = y,
      width = width,
      fill = target,
      alpha = jaccard,
      label = Info
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
      
      panel.spacing = unit(0.05, "lines"),
      
      strip.text.x = element_blank(),
      strip.text.y = element_blank(),
      axis.title.y = element_blank()
    ) +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    scale_fill_manual(values = colors) +
    labs(x = "Position")
  p
  
  
  l <- ggplotly(p, tooltip = c('label') )%>%
    layout(legend = list(orientation = "h", x = 0, y =-0.1))
  l$height <- length(unique(xxx$query)) * nth.best * 8
  #l
  
  htmlwidgets::saveWidget(l, paste0(query, ".html"))
}
