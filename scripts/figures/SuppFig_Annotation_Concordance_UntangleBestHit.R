args <- commandArgs()
path_concordance_by_haplotype_tsv <- args[6]
path_untangle_grounded_tsv <- args[7]
dir_annotation <- args[8]
path_output <- args[9]

x_min <- 0
x_max <- 25500000
widt <- 100
height <- 32

panel_spacing <- 0

nth_best <- 1
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
  'grch38#chr13',
  'HG002#MAT#chr13.prox', 
  'HG002#PAT#chr13.prox',
  'HG002#1#JAHKSE010000070.1',
  'HG002#1#JAHKSE010000214.1',
  'HG002#1#JAHKSE010000225.1',
  'HG002#2#JAHKSD010000014.1',
  'HG002#2#JAHKSD010000065.1',
  'HG002#2#JAHKSD010000105.1',
  'HG002#2#JAHKSD010000106.1',
  'HG002#2#JAHKSD010000126.1',
  'HG002#2#JAHKSD010000186.1'
)
query_to_consider[[length(query_to_consider)+1]] <- c(
  'chm13#chr14',
  'grch38#chr14',
  'HG002#MAT#chr14.prox',
  'HG002#PAT#chr14.prox',
  'HG002#1#JAHKSE010000015.1',
  'HG002#1#JAHKSE010000060.1',
  'HG002#1#JAHKSE010000062.1',
  'HG002#1#JAHKSE010000099.1',
  'HG002#1#JAHKSE010000195.1',
  'HG002#1#JAHKSE010000229.1',
  'HG002#1#JAHKSE010000263.1',
  'HG002#1#JAHKSE010000326.1',
  'HG002#2#JAHKSD010000044.1',
  'HG002#2#JAHKSD010000060.1',
  'HG002#2#JAHKSD010000125.1',
  'HG002#2#JAHKSD010000136.1',
  'HG002#2#JAHKSD010000199.1',
  'HG002#2#JAHKSD010000320.1'
)
query_to_consider[[length(query_to_consider)+1]] <- c(
  'chm13#chr15', 
  'grch38#chr15',
  'HG002#MAT#chr15.prox', 
  'HG002#PAT#chr15.prox',
  'HG002#1#JAHKSE010000018.1',
  'HG002#1#JAHKSE010000073.1',
  'HG002#1#JAHKSE010000097.1',
  'HG002#2#JAHKSD010000064.1',
  'HG002#2#JAHKSD010000094.1',
  'HG002#2#JAHKSD010000123.1'
)
query_to_consider[[length(query_to_consider)+1]] <- c(
  'chm13#chr21', 
  'grch38#chr21',
  'HG002#MAT#chr21.prox', 
  'HG002#PAT#chr21.prox',
  'HG002#1#JAHKSE010000013.1',
  'HG002#1#JAHKSE010000014.1',
  'HG002#2#JAHKSD010000051.1',
  'HG002#2#JAHKSD010000059.1',
  'HG002#2#JAHKSD010000164.1'
)
query_to_consider[[length(query_to_consider)+1]] <- c(
  'chm13#chr22', 
  'grch38#chr22',
  'HG002#MAT#chr22.prox', 
  'HG002#PAT#chr22.prox',
  'HG002#1#JAHKSE010000065.1',
  'HG002#1#JAHKSE010000154.1',
  'HG002#2#JAHKSD010000011.1',
  'HG002#2#JAHKSD010000035.1',
  'HG002#2#JAHKSD010000043.1',
  'HG002#2#JAHKSD010000050.1',
  'HG002#2#JAHKSD010000110.1'
)

c <- read.delim(path_concordance_by_haplotype_tsv)
c$num.different.targets <- as.character(c$num.different.targets)
c[c$num.different.targets != '0' & c$num.different.targets != '1',]$num.different.targets <- '1+'

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
      plot.margin = unit(c(0,8.03,0,1.86), "cm")
    )
  
  
  chr <- paste0('chm13#chr', num_chr)
  
  # Concordance
  c_chr <- c[c$grounded.target == chr,]
  
  p_concordance <- ggplot(
    c_chr,
    aes(
      x = start.pos + (end.pos - start.pos) / 2,
      width = end.pos - start.pos,
      y = haplotype,
      alpha= num.different.targets, fill = grounded.target)
  ) +
    geom_tile() +
    theme_bw() +
    #facet_col(~grounded.target, scales = "free_y", space = "free", labeller = labeller(variable = labels)) +
    theme(
      plot.title = element_text(hjust = 0.5),
      
      text = element_text(size = 32),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.position = "top",
      
      plot.margin = unit(c(0,0.25,0,5.47), "cm"),
      panel.spacing = unit(panel_spacing, "lines")
      #axis.title.y=element_blank()
    ) +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    labs(x = "", y = "", fill="Target") +
    guides(
      fill = "none",
      alpha = guide_legend(title="# different targets", override.aes = list(size=10))
    ) +
    scale_fill_manual(values=colors[[i]])
  
  # Best untangled hit
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
      
      plot.margin = unit(c(0,0,0,0), "cm"),
      
      panel.spacing = unit(panel_spacing, "lines"),
      #panel.border = element_rect(color = "grey", fill = NA, size = 1), #element_blank(),
      
      strip.text.x = element_blank(),
      strip.text.y = element_blank(),
      axis.title.y = element_blank()
    ) +
    scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
    scale_fill_manual(values = colors) +
    labs(x = "Position") +
    guides(
      fill = guide_legend(title="Target", override.aes = list(size=10)),
      alpha = guide_legend(title="Jaccard", override.aes = list(size=10))
    )
  
  
  # Final panel
  p_panel <- ggpubr::ggarrange(
    p_annotation, p_concordance, p_untangle,
    labels=c('', '', ''), font.label = list(size = 40, color = "black", face = "bold", family = NULL),
    heights = c(1.8, 0.7, 2.5),
    legend = "right", # legend position,
    common.legend = F,
    nrow = 3
  )

  ggsave(
    plot = p_panel,
    paste0('SuppFigure', 11 + i, '.pdf'),
    width = width, height = height,
    units = "cm",
    dpi = 100, bg = "white",
    limitsize = FALSE
  )
}
