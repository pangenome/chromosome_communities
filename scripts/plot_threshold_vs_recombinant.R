args <- commandArgs()
path_recombinant_regions_table_sizes_tsv <- args[6]
path_output <- args[7]


library(ggplot2)

x <- read.delim(path_recombinant_regions_table_sizes_tsv, header = F)
x$V3 <- x$V3/1000/1000 # Convert in Mbps
colnames(x) <- c('Max. self coverage', 'Min. estimated identity', 'Mbps')

# Replace 0 with Inf, make it a factor with the specified order
x$`Max. self coverage` <- as.character(x$`Max. self coverage`)
x[x$`Max. self coverage` == 0,]$`Max. self coverage` <- 'Inf'
x$`Max. self coverage` <- factor(x$`Max. self coverage`, levels = c("Inf", "1.1"))
options(scipen=10000) # Disable scientific notation on axes

p <- ggplot(x, aes(x=`Min. estimated identity`, y=Mbps, color=`Max. self coverage`)) +
  geom_line() + 
  theme_bw() + 
  facet_grid(~`Max. self coverage`) +
  scale_x_continuous(breaks = round(seq(min(x$`Min. estimated identity`), max(x$`Min. estimated identity`), by = 0.01), 3)) +
  scale_y_continuous(breaks = round(seq(0, max(x$Mbps), by = 0.5), 3)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(
  plot = p,
  path_output,
  width = 35, height = 15,
  units = "cm",
  dpi = 100, bg = "white",
  limitsize = FALSE
)
