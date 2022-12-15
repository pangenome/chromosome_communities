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
    "CR382285.2", "CR381535.11", "CR381572.5", "CR392039.8","CR381653.11", "CR382287.7", "CR382332.15"
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
  labs(x = "Position (Mbp)", fill="Info", alpha='Estimated identity') +
  scale_alpha(range=c(0.05,1), limits=c(min_id,100)) +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0), breaks=pretty_breaks(n=8)) +
  scale_fill_manual(guide = "none", values = c(
    "chr13-PHRs" = "#F8766D",
    "chr14-PHRs"="#A3A500",
    "chr15-PHRs"="#00BF7D",
    "chr21-PHRs"="#00B0F6",
    "chr22-PHRs"="#E76BF4",
    "chr13-SST1"="#0033FF",
    "chr14-SST1"="#0033FF",
    "chr21-SST1"="#0033FF",
    "chr13-rDNA"="#FF3300",
    "chr14-rDNA"="#FF3300",
    "chr15-rDNA"="#FF3300",
    "chr22-rDNA"="#FF3300",
    "chr21-rDNA"="#FF3300",
    "CR382285.2"="#000000",
    "CR381535.11"="#000000",
    "CR381572.5"="#000000",
    "CR392039.8"="#000000",
    "CR381653.11"="#000000",
    "CR382287.7"="#000000",
    "CR382332.15"="#000000"
    ))
#p
ggsave(
  plot = p,
  file.path(paste0('SupplementaryFigureX9.ROB.clones.id', min_id, '.pdf')),
  width = 45,
  height = 1.5 * (xx %>% pull(info) %>% unique() %>% length()),
  units = "cm",
  dpi = 300, bg = "white",
  limitsize = FALSE
)


library(plotly)
ggplotly(p)
