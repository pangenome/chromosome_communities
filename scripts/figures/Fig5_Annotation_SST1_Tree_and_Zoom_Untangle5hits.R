args <- commandArgs()
path_untangle_grounded_tsv <- args[6]
dir_annotation <- args[7]
path_output <- args[8]

height_bar <- 0.9
panel_spacing <- 0

estimated_identity_threshold <- 0.9
nth_best <- 5
ref_nth_best <- 1


library(ggplot2)
library(ggforce)
library(tidyverse)
library(png)
library(grid)
library(ggpubr)

query_to_consider <- c()
query_to_consider[[length(query_to_consider)+1]] <- c(
  'chm13#chr13',
  'HG002#MAT#chr13.prox',
  'HG002#PAT#chr13.prox',
  'HG01361#2#JAGYYW010000010.1',
  'HG01978#1#JAGYVS010000056.1',
  'HG02486#1#JAGYVM010000043.1',
  'HG03540#2#JAGYVX010000153.1'
)
query_to_consider[[length(query_to_consider)+1]] <- c(
  'chm13#chr14',
  'HG002#MAT#chr14.prox',
  'HG002#PAT#chr14.prox',
  'HG00735#1#JAHBCH010000039.1',
  'HG00741#2#JAHALX010000038.1',
  'HG01106#2#JAHAMB010000013.1',
  'HG01978#1#JAGYVS010000055.1'
)
query_to_consider[[length(query_to_consider)+1]] <- c(
  'chm13#chr21',
  'HG002#MAT#chr21.prox',
  'HG002#PAT#chr21.prox',
  'HG00735#2#JAHBCG010000066.1',
  'HG02886#1#JAHAOU010000106.1',
  'NA18906#1#JAHEOO010000072.1',
  'NA19240#2#JAHEOL010000065.1'
)

x <- read.delim(path_untangle_grounded_tsv) %>%
  filter(
    grounded.target %in% c('chm13#chr13', 'chm13#chr14', 'chm13#chr21') & 
      nth.best <= nth_best & ref.nth.best <= ref_nth_best
    )

# From https://doi.org/10.1093/bioinformatics/btac244
x$estimated_identity <- exp((1.0 + log(2.0 * x$jaccard/(1.0+x$jaccard)))-1.0)

x <- x[x$estimated_identity >= estimated_identity_threshold,]

colors <- c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF4")

chromosomes <- c(13, 14, 21)

# SST1 tree
path_sst1_tree = '/home/guarracino/git/chromosome_communities/data/SST1tree.2022.12.22.png'
img_sst1_tree <- readPNG(path_sst1_tree)
p_sst1_tree <- ggplot() +
  annotation_custom(
    rasterGrob(img_sst1_tree, width = 1, height = 1),
    xmin = - Inf, xmax = Inf,
    ymin = - Inf, ymax = Inf
  ) + theme(
    plot.margin = unit(c(13.5,3,13.5,3), "cm")
  )


delta <- 1000000

p_panels <- c()
for (i in seq_along(chromosomes)){
  num_chr <- chromosomes[[i]]
  
  if (num_chr == 13) {
    x_min <- 12301367 - delta
    x_max <- 12440010 + delta
  } else if (num_chr == 14) {
    x_min <- 6960008 - delta
    x_max <- 6988409 + delta
  } else if (num_chr == 21) {
    x_min <- 9375567 - delta
    x_max <- 9453313 + delta
  }
  
  print(num_chr)

  # Karyotype
  path_karyotype <- file.path(dir_annotation, paste0('genome_browser_chr', num_chr, '_SST1_1Mbps.karyo.png'))
  img_karyo <- readPNG(path_karyotype)
  p_karyotype <- ggplot() +
    annotation_custom(
      rasterGrob(img_karyo, width = 1, height = 1),
      xmin = - Inf, xmax = Inf,
      ymin = - Inf, ymax = Inf
    ) + theme(
      plot.margin = unit(c(0,0.5,0,0.59+7.35), "cm")
    )
  
  # Annotation
  path_annotation <- file.path(dir_annotation, paste0('genome_browser_chr', num_chr, '_SST1_1Mbps_CenSatAnnDense.png'))
  img <- readPNG(path_annotation)
  p_annotation <- ggplot() +
    annotation_custom(
      rasterGrob(img, width = 1, height = 1),
      xmin = - Inf, xmax = Inf,
      ymin = - Inf, ymax = Inf
    ) + theme(
      plot.margin = unit(c(0,0.5,0.1,0.59+3.95), "cm")
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
      x = ref.begin + (ref.end - ref.begin) / 2,
      width = ref.end - ref.begin,
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
      
      text = element_text(size = 24),
      axis.text.x = element_text(size = 13),
      axis.text.y = element_text(size = 13),
      
      legend.title = element_text(size = 30),
      legend.text = element_text(size = 30),
      legend.position = "none",
      
      panel.spacing = unit(panel_spacing, "lines"),
      #panel.border = element_rect(color = "grey", fill = NA, size = 1), #element_blank(),
      
      strip.text.x = element_blank(),
      strip.text.y = element_blank(),
      
      #axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      
      plot.margin = unit(c(0,0.5,0,0), "cm"),
    ) +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    scale_fill_manual(values = colors) +
    labs(x = "Position")+
    guides(
      fill = guide_legend(title="Target", override.aes = list(size=10)),
      alpha = guide_legend(title="Estimated identity", override.aes = list(size=10))
    )

  # Final panel
  p_panel <- ggpubr::ggarrange(
    p_karyotype, p_annotation, p_untangle,
    labels=c('', '', ''), font.label = list(size = 40, color = "black", face = "bold", family = NULL),
    heights = c(0.18, 0.75, 3.8),
    legend = "right", # legend position,
    common.legend = T,
    nrow = 3
  )

  #ggsave(
  #  plot = p_panel,
  #  file.path(path_output, paste0(num_chr, '.png')),
  #  width = width*1, height = height,
  #  units = "cm",
  #  dpi = 100, bg = "white",
  #  limitsize = FALSE
  #)

  p_panels[[length(p_panels)+1]] <- p_panel
}


# Put panels together
p_figure_B <- ggpubr::ggarrange(
  p_panels[[1]], p_panels[[2]], p_panels[[3]],
  labels=c('chr13', 'chr14', 'chr21'),
  font.label = list(size = 35, color = "black", face = "bold", family = NULL),
  hjust=-0.5, vjust=+2.3,
  heights = c(1,1,1),
  
  legend = "right", # legend position,
  common.legend = T,
  nrow = 3, ncol = 1
) + theme(plot.margin = margin(3,0,3,0, "cm")) 

p_figure_AB <- ggpubr::ggarrange(
  p_sst1_tree, p_figure_B,
  labels=c('A', 'B'),
  font.label = list(size = 40, color = "black", face = "bold", family = NULL),
  hjust=-0.1, vjust=+1.1,
  align = c('v'),
  heights = c(0.1,1),
  widths = c(0.4, 0.6),
  legend = "right", # legend position,
  common.legend = T,
  nrow = 1, ncol = 2
)

width <- 90
height <- 70
ggsave(
  plot = p_figure_AB,
  path_output,
  width = width, height = height, # There are 4 NULL plots that have height == height/20
  units = "cm",
  dpi = 200, bg = "white",
  limitsize = FALSE
)

