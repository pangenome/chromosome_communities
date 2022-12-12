args <- commandArgs()
path_fimo_bed <- args[6]
width <- as.numeric(args[7])
path_output <- args[8]

library(ggridges)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(scales) # for pretty_breaks()
library(gtools) # for mixedsort install.packages("gtools")


x <- read.delim(path_fimo_bed, header = F)
colnames(x) <- c('ref', 'ref.begin', 'ref.end', 'label', 'strand', 'qvalue')

x$ref <- gsub('[:].*$', '', x$ref)

x$loq10qvalue <- -log10(x$qvalue)
x[!is.finite(x$loq10qvalue),]$loq10qvalue <- max(x[is.finite(x$loq10qvalue),]$loq10qvalue)


colors <- c("#F8766D", "#A3A500", "#00B0F6")

p <- ggplot(
  x,
  aes(
    x = ref.begin + (ref.end - ref.begin) / 2,
    width = ref.end - ref.begin ,
    y = ordered(label, levels = rev(mixedsort(unique(label)))),
    fill = ref,
    alpha = loq10qvalue
  )
) +
  geom_tile() +
  facet_wrap(~ref, scales = "free", ncol=1, labeller = labeller(variable = labels)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    
    text = element_text(size = 14),
    axis.text.x = element_text(size = 12, angle = 45, vjust = 1.0, hjust=1),
    axis.text.y = element_text(size = 10),
    
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.position = "right",
    
    #panel.spacing = unit(panel_spacing, "lines"),
    #panel.border = element_rect(color = "grey", fill = NA, size = 1), #element_blank(),
    
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    axis.title.y = element_blank(),
    #axis.title.x = element_blank(),
    plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")
  ) +
  ggtitle('PRDM9 motif hits in the SST1 repeat unit') +
  scale_x_continuous(
    #limits = c(x_min, x_max),
    breaks = pretty_breaks(n=12),
    expand = c(0.01, 0.01)) +
  scale_fill_manual(values=colors) +
  labs(x = "Position (bp)", y = '', fill = "Chromosome", alpha = expression("-log"[10]("qvalue")))
#+ scale_alpha_continuous(limits = c(min(1 - x$qvalue), 1))
#p

ggsave(
  plot = p,
  path_output,
  width = width, height = length(unique(x$ref))*7,
  units = "cm",
  dpi = 300, bg = "transparent",
  limitsize = FALSE
)
