args <- commandArgs()
path_untangle_grounded_tsv <- args[6]
dir_annotation <- args[7]
path_output <- args[8]

x_min <- 0
x_max <- 25000000
width <- 90
height <- 36

panel_spacing <- 0.2

nth_best <- 5
ref_nth_best <- 1


library(ggplot2)
library(ggforce)
library(tidyverse)

query_to_consider <- c()
query_to_consider[[length(query_to_consider)+1]] <- c(
  'chm13#chr13', 
  'grch38#chr13',
  'HG002#MAT#chr13.prox', 
  'HG002#PAT#chr13.prox',
  'HG01361#2#JAGYYW010000010.1',
  'HG01978#1#JAGYVS010000056.1',
  'HG02486#1#JAGYVM010000043.1',
  'HG03540#2#JAGYVX010000153.1'
)
query_to_consider[[length(query_to_consider)+1]] <- c(
  'chm13#chr14',
  'grch38#chr14',
  'HG002#MAT#chr14.prox',
  'HG002#PAT#chr14.prox',
  'HG00735#1#JAHBCH010000039.1',
  'HG00741#2#JAHALX010000038.1',
  'HG01978#1#JAGYVS010000055.1',
  'HG02630#1#JAHAOQ010000067.1'
)
query_to_consider[[length(query_to_consider)+1]] <- c(
  'chm13#chr15', 
  'grch38#chr15',
  'HG002#MAT#chr15.prox', 
  'HG002#PAT#chr15.prox',
  'HG00741#2#JAHALX010000004.1',
  'HG02486#2#JAGYVL010000058.1',
  'HG03486#2#JAHEOP010000088.1',
  'NA18906#2#JAHEON010000012.1 '
)
query_to_consider[[length(query_to_consider)+1]] <- c(
  'chm13#chr21', 
  'grch38#chr21',
  'HG002#MAT#chr21.prox', 
  'HG002#PAT#chr21.prox',
  'HG00735#2#JAHBCG010000066.1',
  'HG02886#1#JAHAOU010000106.1',
  'NA18906#1#JAHEOO010000072.1',
  'NA19240#2#JAHEOL010000065.1'
)
query_to_consider[[length(query_to_consider)+1]] <- c(
  'chm13#chr22', 
  'grch38#chr22',
  'HG002#MAT#chr22.prox', 
  'HG002#PAT#chr22.prox',
  'HG00735#1#JAHBCH010000040.1', 
  'HG01361#1#JAGYYX010000045.1', 
  'HG02055#1#JAHEPK010000087.1', 
  'HG03098#1#JAHEPM010000147.1'
)

x <- read.delim(path_untangle_grounded_tsv) %>%
  filter(nth.best <= nth_best & ref.nth.best <= ref_nth_best)

# To have it as numeric column
#x$self.coverage[x$self.coverage == '.'] <- 1
#x$self.coverage <- as.numeric(x$self.coverage)
#x <- x[x$self.coverage <= 1,]

colors <- c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF4")

chromosomes <- c(13, 14, 15, 21, 22)

p_panels <- c()

for (i in seq_along(chromosomes)){
  num_chr <- chromosomes[[i]]
  print(num_chr)
  
  # Annotation
  path_annotation <- file.path(dir_annotation, paste0('hgt_genome_euro_chr', num_chr, '_0_25Mbp.png'))
  img <- readPNG(path_annotation)
  p_annotation <- ggplot() +
    annotation_custom(
      rasterGrob(img, width = 1, height = 1),
      xmin = - Inf, xmax = Inf,
      ymin = - Inf, ymax = Inf
    ) + theme(
      plot.margin = unit(c(0,1,0,0), "cm")
    )
  

  chr <- paste0('chm13#chr', num_chr)
    
  # Apply filters
  xx <- x %>%
    filter(grounded.target == chr & (query %in% query_to_consider[[i]]))

  # Do not consider dedicated annotation bars
  # Do not consider other acros references
  xx <- xx %>%
    filter(!grepl('rDNA', query) & !grepl('centromere', query)) %>%
    filter(!grepl('chr', query) | grepl(paste0('chr', num_chr), query))
  
  # To group by query
  xx$query.hacked <- paste(xx$query, xx$nth.best, sep = "-")
  
  xx <- xx %>%
    arrange(query.hacked)
  
  p_untangle <- ggplot(
    xx,
    aes(
      x = ref.begin + (ref.end - ref.begin) / 2, width = ref.end - ref.begin - 200,
      y = ordered(query.hacked, levels = rev(unique(query.hacked))),
      fill = target,
      alpha = jaccard
    )
  ) +
    geom_tile() +
    #ggtitle(paste(chr, title)) +
    facet_grid(query~grounded.target, scales = "free_y", space = "free", labeller = labeller(variable = labels)) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      
      text = element_text(size = 32),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.position = "top",
      
      panel.spacing = unit(panel_spacing, "lines"),
      #panel.border = element_rect(color = "grey", fill = NA, size = 1), #element_blank(),
      
      strip.text.x = element_blank(),
      strip.text.y = element_blank(),
      axis.title.y = element_blank()
    ) +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    scale_fill_manual(values = colors) +
    labs(x = "Position")+
    guides(
      fill = guide_legend(title="Target", override.aes = list(size=10)),
      alpha = guide_legend(title="Jaccard", override.aes = list(size=10))
    )

  # Final panel
  p_panel <- ggpubr::ggarrange(
    p_annotation, p_untangle,
    labels=c('', '', ''), font.label = list(size = 40, color = "black", face = "bold", family = NULL),
    heights = c(1.8, 4.4),
    legend = "right", # legend position,
    common.legend = T,
    nrow = 2
  )

  p_panels[[length(p_panels)+1]] <- p_panel
}


# Put panels together
p_figure <- ggpubr::ggarrange(
  p_panels[[1]], p_panels[[2]], p_panels[[3]], p_panels[[4]], p_panels[[5]],
  labels=c('A', 'B', 'C', 'D', 'E'), font.label = list(size = 40, color = "black", face = "bold", family = NULL),
  heights = c(1,1,1,1,1),
  legend = "right", # legend position,
  common.legend = T,
  nrow = 5, ncol = 1
)

ggsave(
  plot = p_figure,
  file.path(path_output, 'Figure5.pdf'),
  width = width*1, height = height*5,
  units = "cm",
  dpi = 100, bg = "white",
  limitsize = FALSE
)
ggsave(
  plot = p_figure,
  file.path(path_output, 'Figure5.png'),
  width = width*1, height = height*5,
  units = "cm",
  dpi = 100, bg = "white",
  limitsize = FALSE
)
