library(ggplot2)
library(ggforce)
library(tidyverse)
library(scales) # for pretty_breaks()


options(scipen = 9)
x <- read.delim(
  '/home/guarracino/Downloads/Pangenomics/chromosome_communities/Review1/robertsonian_translocation/chrACRO_7-Dec-22_PHRs+PMC4257996+SST1+rDNA.bed',
  header = F
)
colnames(x) <- c('chrom', 'begin', 'end', 'identity', 'info')

x$identity <- as.numeric(x$identity)

x[x$chrom == 'chm13#chr13' & x$info == 'SST1',]$info <- 'chr13-SST1'
x[x$chrom == 'chm13#chr14' & x$info == 'SST1',]$info <- 'chr14-SST1'
x[x$chrom == 'chm13#chr21' & x$info == 'SST1',]$info <- 'chr21-SST1'

x[x$chrom == 'chm13#chr13' & x$info == 'rDNA',]$info <- 'chr13-rDNA'
x[x$chrom == 'chm13#chr14' & x$info == 'rDNA',]$info <- 'chr14-rDNA'
x[x$chrom == 'chm13#chr15' & x$info == 'rDNA',]$info <- 'chr15-rDNA'
x[x$chrom == 'chm13#chr21' & x$info == 'rDNA',]$info <- 'chr21-rDNA'
x[x$chrom == 'chm13#chr22' & x$info == 'rDNA',]$info <- 'chr22-rDNA'

x[x$chrom == 'chm13#chr13' & x$info == 'PHRs',]$info <- 'chr13-PHRs'
x[x$chrom == 'chm13#chr14' & x$info == 'PHRs',]$info <- 'chr14-PHRs'
x[x$chrom == 'chm13#chr15' & x$info == 'PHRs',]$info <- 'chr15-PHRs'
x[x$chrom == 'chm13#chr21' & x$info == 'PHRs',]$info <- 'chr21-PHRs'
x[x$chrom == 'chm13#chr22' & x$info == 'PHRs',]$info <- 'chr22-PHRs'




x$info <- factor(
  x$info,
  levels = rev(c(
    "chr13-PHRs", "chr14-PHRs", "chr15-PHRs", "chr21-PHRs", "chr22-PHRs",
    "chr13-SST1", "chr14-SST1", "chr21-SST1",
    "chr13-rDNA", "chr14-rDNA", "chr15-rDNA", "chr21-rDNA", "chr22-rDNA",
    "CR382285.2", "CR381535.11", "CR381572.5", "CR392039.8","CR381653.11", "CR382287.7", "CR382332.15", "CR381670.5"
    ))
)


x_min <- 0
x_max <- 18000000#max(x %>% filter(chrom %in% c('chm13#chr13', 'chm13#chr14', 'chm13#chr15', 'chm13#chr21', 'chm13#chr22') & identity >= 96) %>% pull(end))
min_id <- 90 # 90 and 99 min(x$identity)

xx <- x %>% filter(chrom %in% c('chm13#chr13', 'chm13#chr14', 'chm13#chr15', 'chm13#chr21', 'chm13#chr22') & identity >= min_id)
p <- ggplot(
  xx,
  aes(
    x = begin + (end - begin) / 2,
    width = end - begin,
    #y = ordered(info, levels = rev(unique(info))),
    y = info,
    fill = info,
    alpha = identity
  )
) +
  geom_tile() +
  #ggtitle(paste(chr, title)) +
  #facet_grid(~chrom, scales = "free_y", space = "free", labeller = labeller(variable = labels)) +
  #facet_grid( ~chrom, ncol = 1, scales = "free_y", space = "free_y") +
  facet_grid(rows = vars(chrom), scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    
    text = element_text(size = 18),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.position = "top",
    
    #panel.spacing = unit(panel_spacing, "lines"),
    #panel.border = element_rect(color = "grey", fill = NA, size = 1), #element_blank(),
    
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.y = element_blank(),
    
    plot.margin = unit(c(0.1,1.3,0.1,0.3), "cm"),
  ) +
  labs(x = "Position (bp)", fill="Info", alpha='Estimated identity') +
  scale_alpha(range=c(0.05,1), limits=c(min_id,100)) +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0), breaks=pretty_breaks(n=8)) +
  scale_fill_manual(guide = "none", values = c(
    "chr13-PHRs" = "#F8766D",
    "chr14-PHRs"="#A3A500",
    "chr15-PHRs"="#00BF7D",
    "chr21-PHRs"="#00B0F6",
    "chr22-PHRs"="#E76BF4",
    "chr13-SST1"="#000000",
    "chr14-SST1"="#000000",
    "chr21-SST1"="#000000",
    "chr13-rDNA"="#000000",
    "chr14-rDNA"="#000000",
    "chr15-rDNA"="#000000",
    "chr22-rDNA"="#000000",
    "chr21-rDNA"="#000000",
    "CR382285.2"="#00FF00",
    "CR381535.11"="#FF0000",
    "CR381572.5"="#00FF00",
    "CR392039.8"="#FF0000",
    "CR381653.11"="#FF0000",
    "CR382287.7"="#00FF00",
    "CR382332.15"="#FF0000",
    "CR381670.5"="#FF0000"
    ))
#p
ggsave(
  plot = p,
  file.path(paste0('SupplementaryFigure28.ROB.clones.id', min_id, '.pdf')),
  width = 45,
  height = 1.5 * (xx %>% pull(info) %>% unique() %>% length()),
  units = "cm",
  dpi = 300, bg = "white",
  limitsize = FALSE
)


# library(plotly)
#ggplotly(p)



# Panel chr14 reversed and chr21
library(ggpubr)

delta <- 1000000

p_panels <- c()

for (num_chr in c(14, 21)) {
  print(num_chr)
  z <- xx %>% filter(chrom == paste0('chm13#chr', num_chr) & identity >= min_id)
  
  if (num_chr == 14) {
    z$info <- factor(
      z$info,
      levels = rev(c(
        "CR382285.2", "CR381535.11", "CR381572.5", "CR392039.8","CR381653.11", "CR382287.7", "CR382332.15", "CR381670.5",
        "chr13-rDNA", "chr14-rDNA", "chr15-rDNA", "chr21-rDNA", "chr22-rDNA",
        "chr13-SST1", "chr14-SST1", "chr21-SST1",
        "chr13-PHRs", "chr14-PHRs", "chr15-PHRs", "chr21-PHRs", "chr22-PHRs"
      ))
    )
  }
  
  
  p <- ggplot(
    z,
    aes(
      x = begin + (end - begin) / 2,
      width = end - begin,
      #y = ordered(info, levels = rev(unique(info))),
      y = info,
      fill = info,
      alpha = identity
    )
  ) +
    geom_tile() +
    #ggtitle(paste(chr, title)) +
    #facet_grid(~chrom, scales = "free_y", space = "free", labeller = labeller(variable = labels)) +
    #facet_grid( ~chrom, ncol = 1, scales = "free_y", space = "free_y") +
    facet_grid(rows = vars(chrom), scales = "free_y", space = "free_y") +
    theme_bw() +
    labs(x = "Position (bp)", fill="Info", alpha='Estimated identity') +
    scale_alpha(range=c(0.05,1), limits=c(min_id,100)) +
    scale_fill_manual(guide = "none", values = c(
      "chr13-PHRs" = "#F8766D",
      "chr14-PHRs"="#A3A500",
      "chr15-PHRs"="#00BF7D",
      "chr21-PHRs"="#00B0F6",
      "chr22-PHRs"="#E76BF4",
      "chr13-SST1"="#000000",
      "chr14-SST1"="#000000",
      "chr21-SST1"="#000000",
      "chr13-rDNA"="#000000",
      "chr14-rDNA"="#000000",
      "chr15-rDNA"="#000000",
      "chr22-rDNA"="#000000",
      "chr21-rDNA"="#000000",
      "CR382285.2"="#00FF00",
      "CR381535.11"="#FF0000",
      "CR381572.5"="#00FF00",
      "CR392039.8"="#FF0000",
      "CR381653.11"="#FF0000",
      "CR382287.7"="#00FF00",
      "CR382332.15"="#FF0000",
      "CR381670.5"="#FF0000"
    ))
  if (num_chr == 14) {
    p <- p + scale_x_reverse(limits = c(x_max - 1610000, x_min), expand = c(0, 0), breaks=pretty_breaks(n=8))
    p <- p +
      theme(
        plot.title = element_text(hjust = 0.5),
        
        text = element_text(size = 21),
        axis.text.x =  element_text(size = 21),
        axis.text.y = element_text(size = 19),
        
        legend.title = element_text(size = 21),
        legend.text = element_text(size = 21),
        legend.position = "top",
        
        #panel.spacing = unit(panel_spacing, "lines"),
        #panel.border = element_rect(color = "grey", fill = NA, size = 1), #element_blank(),
        
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        
        plot.margin = unit(c(0.1,0.3,0.1,0.3), "cm"),
        #plot.margin = unit(c(0.1,6.45,0.1,0.3), "cm"),
      ) +
      # Centromere
      geom_vline(xintercept=12776582, color='black', alpha=0.3, size=1) +
      geom_vline(xintercept=10067897, color='black', alpha=0.3, size=1) +
      # SST1 +- 1Mbp
      geom_vline(xintercept=6960008 - delta, color='black', alpha=1, size=1) +
      geom_vline(xintercept=6988409 + delta, color='black', alpha=1, size=1)
      #annotate("rect",
      #            xmin = 12776582, xmax = 10067897,
      #            ymin = 1 - 0.6,
      #            ymax = 11 + 0.6,
      #            fill = "#f8ec32", alpha = .2, color = "#f8ec32", size = 0)
  } else {
    p <- p + scale_x_continuous(limits = c(x_min, x_max - 1610000), expand = c(0, 0), breaks=pretty_breaks(n=8))
    p <- p +
      theme(
        plot.title = element_text(hjust = 0.5),
        
        text = element_text(size = 21),
        axis.text.x = element_text(size = 21),
        axis.text.y = element_text(size = 19),
        
        legend.title = element_text(size = 21),
        legend.text = element_text(size = 21),
        legend.position = "top",
        
        #panel.spacing = unit(panel_spacing, "lines"),
        #panel.border = element_rect(color = "grey", fill = NA, size = 1), #element_blank(),
        
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        axis.title.y = element_blank(),
        
        plot.margin = unit(c(0.1,0.3,0.1,0.3), "cm"),
      ) +
      # Centromere
      geom_vline(xintercept=10816799, color='black', alpha=0.3, size=1) +
      geom_vline(xintercept=11340698, color='black', alpha=0.3, size=1) +
      # SST1 +- 1Mbp
      geom_vline(xintercept=9375567 - delta, color='black', alpha=1, size=1) +
      geom_vline(xintercept=9453313 + delta, color='black', alpha=1, size=1)
      #annotate("rect",
      #               xmin = 10816799, xmax = 11340698,
      #               ymin = 1 - 0.6,
      #               ymax = 10 + 0.6,
      #               fill = "#f8ec32", alpha = .2, color = "#f8ec32", size = 0)
  }
  
  p_panels[[length(p_panels)+1]] <- p
}

# Approximate centromeres
#chm13#chr13	13941594	17573031	chr13_q#F8766D
#chm13#chr14	10067897	12776582	chr14_q#A3A500
#chm13#chr15	15412039	17709803	chr15_q#00BF7D
#chm13#chr21	10816799	11340698	chr21_q#00B0F6
#chm13#chr22	12784333	16188127	chr22_q#E76BF4

p_panel <- ggpubr::ggarrange(
  p_panels[[1]], p_panels[[2]],
  labels=c('', ''), font.label = list(size = 40, color = "black", face = "bold", family = NULL),
  heights = c(1, 1),
  legend = "none", # legend position,
  common.legend = T,
  nrow = 2
)
p_panel

ggsave(
  plot = p_panel,
  file.path(paste0('Figure4C.ROB.clones.id', min_id, '.chr14inv_and_chr21.pdf')),
  width = 45,
  height = 0.75 * (xx %>% pull(info) %>% unique() %>% length()),
  units = "cm",
  dpi = 300, bg = "white",
  limitsize = FALSE
)

