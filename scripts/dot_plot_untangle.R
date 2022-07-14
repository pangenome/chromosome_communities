#f=chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#ACRO.e50000.m1000.j0.n100.fixed.bed.gz
#(zcat $f | head -n 1; \
#  zgrep '^HG002#MAT\|^HG002#PAT' $f | awk '$8 == "-" { x=$6; $6=$5; $5=x; } { print }' | awk '$10 == 1') | pigz -c > HG002.verkko.bed.gz

path_untangle_bed <- '/home/guarracino/HG002.verkko.bed.gz'

library(ggplot2)
#library(knitr)
x <- read.delim(path_untangle_bed, sep = '\t', header = T, comment.char = "$")
colnames(x)[1] <- 'query.name'

x$query.name <- gsub(":.*", "", x$query.name)
x$query.name <- gsub("#J.*", "", x$query.name)

x <- x[x$ref.start < 25000000 & x$ref.end <=25000000,]

x_subset <- subset(x, query.name %in% c(
  "HG002#MAT#chr13.prox"
)
)

p <- ggplot(
  x_subset, aes(x = query.start, xend = query.end, y = ref.start, yend = ref.end, color=ref.name)) +
  geom_segment(size = 1.5) +
  facet_grid(query.name ~ ref.name) +
  coord_fixed() +
  theme_bw() +
  theme(
    text = element_text(size = 12.6),
    axis.text.x = element_text(size = 12, angle = 90),
    axis.text.y = element_text(size = 12),
    
    #panel.grid.minor = element_line(size = 0.125),
    #panel.grid.major = element_line(size = 0.25)
  ) +
  xlab("Query start") +
  ylab("Reference start")  + theme(aspect.ratio=1)
p

# C4A annotation
# https://stackoverflow.com/questions/65013846/why-geom-rect-looking-different-across-facets
p <- p +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 31982057 - 31972046, ymax = 32002681 - 31972046, fill = "#EE2222", alpha = .2, color = "#444444", size = 0.1)
#https://stackoverflow.com/questions/11889625/annotating-text-on-individual-facet-in-ggplot2
C4A_ann_text <- data.frame(
  `query.start` = 5000, `query.end` = 0, `ref.start` = 20000, `ref.end` = 0, `query.name` = factor("grch38#chr6", levels = unique(x_subset$query.name))
)
p <- p + geom_text(data = C4A_ann_text, label = 'C4A', size = 4)

# C4B annotation
# https://stackoverflow.com/questions/65013846/why-geom-rect-looking-different-across-facets
p <- p +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 32014795 - 31972046, ymax = 32035418 - 31972046, fill = "#2222EE", alpha = .2, color = "#444444", size = 0.1)
C4B_ann_text <- data.frame(
  `query.start` = 5000, `query.end` = 0, `ref.start` = 52500, `ref.end` = 0, `query.name` = factor("grch38#chr6", levels = unique(x_subset$query.name))
)
p <- p + geom_text(data = C4B_ann_text, label = 'C4B', size = 4)


filename <- paste0(gsub("\\.", "_", path_untangle_bed), '.pdf')
#knitr::plot_crop(filename)
ggsave(plot = p, filename, width = 32, height = 8, units = "cm", dpi = 300, bg = "transparent")