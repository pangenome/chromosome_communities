library(karyoploteR)

custom.genome <- toGRanges("/home/guarracino/git/chromosome_communities/data/chm13.acros.tsv")

#https://github.com/marbl/CHM13/issues/47
#https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/CHM13/assemblies/annotation/
#Add header and set the extension to `txt` to avoid problems with `toGRanges`
custom.cytobands <- toGRanges("/home/guarracino/git/chromosome_communities/data/chm13v2.0_cytobands_allchrs.modified.txt")

#Open the device
pdf("Acros.cytobands.pq_arms.pdf", width = 1600, height = 700, bg = "transparent")

#Draw
p_arms <- toGRanges(
  data.frame(
    chr=c("chr13", "chr14", "chr15", "chr21", "chr22"),
    start=c(1, 1, 1, 1, 1),
    end=c(12941594, 9067897, 14412039, 9816799, 11784333)
  )
)
q_arms <- toGRanges(
  data.frame(
    chr=c("chr13", "chr14", "chr15", "chr21", "chr22"),
    start=c(18573031, 13776582, 18709803, 12340698, 17188127),
    end=c(113566686, 101161492, 99753195, 45090682, 51324926)
  )
)
plot.params <- getDefaultPlotParams(plot.type = 2)
plot.params$data1height <- 100

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, plot.params = plot.params, cex=2)
kpPlotRegions(kp, p_arms, col="#CCFFAA", r1=0.3)
kpPlotRegions(kp, q_arms, col="#CCFFAA", r1=0.3)
kpAddBaseNumbers(kp, tick.dist=5000000, cex=1)

#Close the device
dev.off()


gains <- toGRanges(
  data.frame(
    chr=c("chr13", "chr14", "chr15", "chr21", "chr22"),
    start=c(1, 1, 1, 1, 1),
    end=c(25000000, 25000000, 25000000, 25000000, 25000000)
  )
)

kp <- plotKaryotype(genome = custom.genome, cytobands = custom.cytobands, cex=2)
kpPlotRegions(kp, gains, col="#CCCCCC", r1=0.3)
kpAddBaseNumbers(kp, tick.dist=5000000, cex=1)


if (FALSE) {
  kp <- plotKaryotype(genome="hg38", plot.type=1, chromosomes=c("chr13", "chr14", "chr15", "chr21", "chr22"))
  kpPlotRegions(kp, gains, col="#CCFFAA")
}