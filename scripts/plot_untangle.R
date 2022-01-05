args <- commandArgs()
path_untangle_grounded_tsv <- args[6]

library(ggplot2)
x <- read.delim(path_untangle_grounded_tsv)
p <- ggplot(x, aes(x = ref.begin + (ref.end - ref.begin) / 2, width = ref.end - ref.begin - 200, y = query, fill = target)) +
  geom_tile()
ggsave(plot = p, paste0(path_untangle_grounded_tsv, '.pdf'), width = 40, height = 65, units = "cm", dpi = 300, bg = "transparent")
