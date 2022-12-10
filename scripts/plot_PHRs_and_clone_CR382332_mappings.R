require(ggplot2)
require(gggenes)
x <- read.delim('~/SST1_region.gggenes.tsv')
x <- x[x$molecule %in% c('chm13#chr13:11501367-13240010', 'chm13#chr14:6160008-7788409', 'chm13#chr21:8575567-10253313'),]
p <- ggplot(x, aes(xmin=start, xmax=end, y=molecule, fill=gene, forward=strand)) + geom_gene_arrow()
ggsave(plot=p, '~/SST1_region.gggenes.subset.png', height=1.5, width=15)


library(ggplot2)
library(ggforce)
library(tidyverse)
library(plotly)
library(scales) # for pretty_breaks()


options(scipen = 9)
x <- read.delim(
  '/home/guarracino/Downloads/Pangenomics/chromosome_communities/Review1/robertsonian_translocation/chrACRO_7-Dec-22_PHRs+CR382332.bed',
  header = F
)
colnames(x) <- c('chrom', 'begin', 'end', 'identity', 'info')

x$identity <- as.numeric(x$identity)


x[x$chrom == 'chm13#chr13' & x$info != 'CR382332',]$info <- 'chr13-PHRs'
x[x$chrom == 'chm13#chr14' & x$info != 'CR382332',]$info <- 'chr14-PHRs'
x[x$chrom == 'chm13#chr15' & x$info != 'CR382332',]$info <- 'chr15-PHRs'
x[x$chrom == 'chm13#chr21' & x$info != 'CR382332',]$info <- 'chr21-PHRs'
x[x$chrom == 'chm13#chr22' & x$info != 'CR382332',]$info <- 'chr22-PHRs'

x_min <- 0
x_max <- 25000000

p <- ggplot(
  x %>% filter(chrom %in% c('chm13#chr13', 'chm13#chr14', 'chm13#chr15', 'chm13#chr21', 'chm13#chr22')),
  aes(
    x = begin + (end - begin) / 2, width = end - begin,
    y = ordered(info, levels = rev(unique(info))),
    fill = info,
    alpha = identity
  )
) +
  geom_tile() +
  #ggtitle(paste(chr, title)) +
  #facet_grid(~chrom, scales = "free_y", space = "free", labeller = labeller(variable = labels)) +
  facet_wrap( ~chrom, ncol = 1, scales = "free_y") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    
    text = element_text(size = 20),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.position = "none",
    
    #panel.spacing = unit(panel_spacing, "lines"),
    #panel.border = element_rect(color = "grey", fill = NA, size = 1), #element_blank(),
    
    #strip.text.x = element_blank(),
    #strip.text.y = element_blank(),
    axis.title.y = element_blank(),
    
    plot.margin = unit(c(0.1,1.3,0.1,0.3), "cm"),
  ) +
  labs(x = "Position", fill="Info") +
  scale_alpha(range=c(0.05,1), limits=c(96,100)) +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0), breaks=pretty_breaks(n=8)) +
  scale_fill_manual(values = c(
    "chr13-PHRs" = "#F8766D",
    "chr14-PHRs"="#A3A500",
    "chr15-PHRs"="#00BF7D",
    "chr21-PHRs"="#00B0F6",
    "chr22-PHRs"="#E76BF4",
    "CR382332"="#000000"
    ))
#p
ggsave(
  plot = p,
  file.path(paste0('SupplementaryFigureX9.ROB.pdf')),
  width = 30, height = 13,
  units = "cm",
  dpi = 300, bg = "white",
  limitsize = FALSE
)



ggplotly(p)
