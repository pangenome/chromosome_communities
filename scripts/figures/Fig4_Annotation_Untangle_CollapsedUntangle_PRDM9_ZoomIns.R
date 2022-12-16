args <- commandArgs()
path_support_dedup_tsv <- args[6]
path_entropy_tsv <- args[7]
dir_annotation <- args[8]
path_output <- args[9]

x_min <- 0
x_max <- 25000000
width <- 70
height <- 40


# library
library(ggridges)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(scales) # for pretty_breaks()

library(png)
library(grid)
library(ggpubr)

options(scipen = 9)


num_chr <- 13
estimated_identity_threshold <- 0.9
height_bar <- 0.7
panel_spacing <- 0

# Karyotype
path_karyotype <- file.path(dir_annotation, paste0('hgt_genome_euro_chr', num_chr, '_0_25Mbp.karyo.png'))
img_karyo <- readPNG(path_karyotype)
p_karyotype <- ggplot() +
  annotation_custom(
    rasterGrob(img_karyo, width = 1, height = 1),
    xmin = - Inf, xmax = Inf,
    ymin = - Inf, ymax = Inf
  ) + theme(
    plot.margin = unit(c(0,0.5,0,7.32), "cm")
  )


# Annotation
path_annotation <- file.path(dir_annotation, paste0('hgt_genome_euro_chr', num_chr, '_0_25Mbp.png'))
img_anno <- readPNG(path_annotation)
p_annotation <- ggplot() +
  annotation_custom(
    rasterGrob(img_anno, width = 1, height = 1),
    xmin = - Inf, xmax = Inf,
    ymin = - Inf, ymax = Inf
  ) + theme(
    plot.margin = unit(c(0,0.5,0.1,2.12), "cm")
  )




chr <- paste0('chm13#chr', num_chr)

query_to_consider <- read.delim(path_query_to_consider, header = F)

u <- read.delim(path_untangle_grounded_tsv) %>%
  filter(
    query %in% query_to_consider$V1 & grounded.target == chr & nth.best <= 1 & ref.nth.best <= 1
    )
rm(query_to_consider)

# From https://doi.org/10.1093/bioinformatics/btac244
u$estimated_identity <- exp((1.0 + log(2.0 * u$jaccard/(1.0+u$jaccard)))-1.0)

u <- u[u$estimated_identity >= estimated_identity_threshold,]



# Do not consider dedicated annotation bars
# Do not consider other acros references
u <- u %>%
  filter(!grepl('rDNA', query) & !grepl('centromere', query)) %>%
  filter(!grepl('chr', query) | grepl(paste0('chr', num_chr), query))

u <- u %>%
  arrange(query)


p_untangle <- ggplot(
  u,
  aes(
    x = ref.begin + (ref.end - ref.begin) / 2, width = ref.end - ref.begin,
    y = ordered(query, levels = rev(unique(query))),
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
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 13),
    
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 30),
    legend.position = "top",
    
    panel.spacing = unit(panel_spacing, "lines"),
    #panel.border = element_rect(color = "grey", fill = NA, size = 1), #element_blank(),
    
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    
    plot.margin = unit(c(0,0.5,0,0), "cm"),
  ) + scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = pretty_breaks(n=20),
    expand = c(0.0, 0.0)
  ) +
  scale_fill_manual(values = c("#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF4")) +
  labs(x = "", fill="Target", alpha="Estimated identity") 
#+ scale_alpha_discrete(range = c(0.3, 1))# + scale_x_reverse()
#ggsave(plot = p_untangle, paste0(path_untangle_grounded_all_tsv, '.pdf'), width = width, height = height, units = "cm", dpi = 100, bg = "transparent", limitsize = FALSE)






# Collapsed untangled output
s <- read.delim(path_support_dedup_tsv)
colnames(s) <- c("ground.target", "start", "end", "chr13", "chr14", "chr15", "chr21", "chr22")


# Entropy
e <- read.delim(path_entropy_tsv)
e[e$shannon_div_index == -1, ]$shannon_div_index <- NA

# Compute average SDI by window by ground.target
e_average <- e %>% 
  dplyr::group_by(ground.target, start.pos, end.pos) %>% 
  dplyr::summarize(average_sdi = mean(shannon_div_index, na.rm = TRUE))



# Apply filters
s_chr <- s[s$ground.target == chr, ]
s_chr_long <- pivot_longer(s_chr, chr13:chr22, "chromosome")
s_chr_long$chromosome <- paste0('chm13#', s_chr_long$chromosome)

if (num_chr == 13) {
  colors <- c("#F8766D")
} else if (num_chr == 14) {
  colors <- c("#A3A500")
} else if (num_chr == 15) {
  colors <- c("#00BF7D")
} else if (num_chr == 21) {
  colors <- c("#00B0F6")
} else if (num_chr == 22) {
  colors <- c("#E76BF4")
}else {
  colors <- c("#000000")
}

p_collapsed_untangle <- ggplot(s_chr_long, aes(x=start, y=value, color=chromosome)) +
  geom_step(alpha=0.5) +
  scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = pretty_breaks(n=20),
    expand = c(0.0, 0.0)
  ) +
  scale_y_continuous(
    limits = c(-max(10, max(s_chr_long$value))*0.05, max(10, max(s_chr_long$value))),
    breaks=pretty_breaks(n=6)
    ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    
    text = element_text(size = 24),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.position = "top",
    
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    
    axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    
    plot.margin = unit(c(0,0.5,0,4.96), "cm")
  ) + labs(
    x = paste(''),
    y = paste('Contigs\n')
  ) +
  guides(
    colour = guide_legend(title="Chromosome", override.aes = list(size=10))
  )
# Gray-out regions with missing aligned contigs
p_collapsed_untangle <- p_collapsed_untangle +
  annotate("rect",
           xmin = 0, xmax = s_chr[1,]$end,
           ymin = 0, ymax = max(10, max(s_chr_long$value)),
           fill = "#444444", alpha = .1, color = "#ffffff", size = 0.1)
# xxxxxxxxxxxxxx
p_collapsed_untangle <- p_collapsed_untangle +
  annotate("rect",
           xmin = 16300000, xmax = 17700000,
           ymin = -max(10, max(s_chr_long$value))*0.05, ymax = max(10, max(s_chr_long$value)),
           fill = "#dc6539", alpha = .0, color = "#ff0000", size = 0.8)
# xxxxxxxxxxxxxx
p_collapsed_untangle <- p_collapsed_untangle +
  annotate("rect",
           xmin = 10700000, xmax = 14100000,
           ymin = -max(10, max(s_chr_long$value))*0.05, ymax = max(10, max(s_chr_long$value)),
           fill = "#dc6539", alpha = .0, color = "#ff0000", size = 0.8)

# Apply filters
e_chr_average <- e_average[e_average$ground.target == chr,]

p_entropy_average <- ggplot(e_chr_average, aes(x=start.pos, y=average_sdi, color=ground.target)) +
  geom_step(alpha=0.8) +
  scale_x_continuous(
    limits = c(x_min, x_max),
    breaks = pretty_breaks(n=20),
    expand = c(0.0, 0.0)
  ) +
  scale_y_continuous(
    limits = c(-min(2, max(e_chr_average$average_sdi, na.rm = T))*0.05, min(2, max(e_chr_average$average_sdi, na.rm = T))),
    breaks=pretty_breaks(n=6)
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    
    text = element_text(size = 24),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14),
    
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.position = "top",
    
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    
    axis.title.x = element_blank(),
    #axis.title.y = element_blank(),
    
    plot.margin = unit(c(0,0.5,0,4.825), "cm")
  ) + labs(
    x = paste('Position'),
    y = paste('Entropy\n')
  ) +
  scale_color_manual(values=colors) +
  guides(colour = guide_legend(override.aes = list(size=10)))
# Gray-out regions with missing aligned contigs
p_entropy_average <- p_entropy_average +
  annotate("rect",
           xmin = 0, xmax = e_chr_average[!is.na(e_chr_average$average_sdi),][1,]$start.pos,
           ymin = 0, ymax = min(2, max(e_chr_average$average_sdi, na.rm = T)),
           fill = "#444444", alpha = .1, color = "#ffffff", size = 0.1)



x <- read.delim(path_fimo_window_bed, header = F)
colnames(x) <- c('chrom', 'ref.begin', 'ref.end', 'hits')

#x$chrom <- gsub('[:].*$','', x$chrom)

colors <- c("#F8766D")

xx <- x[x$chrom == chr & x$ref.begin >= x_min & x$ref.end <= x_max,]
p_PRDM9 <- ggplot(xx, aes(
  x = (ref.begin + (ref.end - ref.begin) / 2) / 1000000, width = ref.end - ref.begin,
  y=hits,
  color=chrom)
  ) +
  geom_step(size=0.3) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    
    text = element_text(size = 24),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.position = "top",
    
    strip.text.x = element_blank(),
    strip.text.y = element_blank(),
    plot.margin = unit(c(0,0.5,0,4.685), "cm")
  ) + scale_x_continuous(
    limits = c(x_min / 1000000, x_max / 1000000),
    breaks = pretty_breaks(n=20),
    expand = c(0.0, 0.0)
    ) + labs(
      x = paste('Position (Mbp)'),
      y = paste('PRDM9\nmotif hits'),
      color = "Chromosome"
    ) + scale_y_continuous(
      #limits = c(0, max(10, max(xx$hits))),
      breaks=pretty_breaks(n=6)
    ) +
  scale_color_manual(values=colors) +
  guides(colour = guide_legend(override.aes = list(size=10)))

xx_min <- 12301367 - 150000
xx_max <- 12440010 + 150000
xxx_max <- max(xx[xx$ref.begin >= xx_min & xx$ref.end <= xx_max,]$hits)
p_PRDM9 <- p_PRDM9 +
  annotate("rect",
           xmin = xx_min / 1000000, xmax = xx_max / 1000000,
           ymin = -xxx_max*0.1, ymax = xxx_max*1.2,
           fill = "#dc6539", alpha = .0, color = "#ff0000", size = 0.8)

# Final panel
p_panel <- ggpubr::ggarrange(
  p_karyotype, p_annotation, p_untangle, p_collapsed_untangle, p_entropy_average, p_PRDM9,
  labels=c('', 'A', '',  'B', 'C', 'D'), font.label = list(size = 40, color = "black", face = "bold", family = NULL),
  heights = c(0.6, 2.5, 7.5, 2.2, 2.2, 2.9),
  legend = "right", # legend position,
  common.legend = T,
  nrow = 6
)
#font.label = list(size = 50, color = "black", face = "bold", family = NULL), hjust=-0.1, vjust=+2,
ggsave(
  plot = p_panel,
  'fig4.png',
  width = width, height = height,
  units = "cm",
  dpi = 100, bg = "white",
  limitsize = FALSE
)
#ggsave(plot = p, path_output, width = width, height = length(unique(x$chrom))*4, units = "cm", dpi = 300, bg = "transparent", limitsize = FALSE)




