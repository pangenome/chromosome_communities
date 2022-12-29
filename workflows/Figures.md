# Figures

## Main figures

### Figure 4

#### Revision 1

```shell
path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.ALL.e50000.m1000.grounded.pq_touching.reliable.tsv.gz
Rscript /lizardfs/guarracino/chromosome_communities/scripts/figures/Fig4_Annotation_Untangle_CollapsedUntangle_PRDM9_ZoomIns.R \
  $path_grounded_pq_touching_reliable_tsv_gz \
  <(zgrep '^HG002#1\|^HG002#2' -v $path_grounded_pq_touching_reliable_tsv_gz | sed '1d' | cut -f 1 | sort | uniq) \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.support.dedup.eid0900.n1.nref1.tsv.gz \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.entropy.by_contig.eid0900.w50000.n1.nref1.tsv \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.entropy_match_order.eid0900.n5.tsv \
  /lizardfs/guarracino/chromosome_communities/recombination_hotspots/chm13v2.PRDM9.w20000.bed
  /lizardfs/guarracino/chromosome_communities/data/annotation/ \
  13 \
  ~/Figure5.pdf
```


path_query_to_consider <- '/home/guarracino/Downloads/Pangenomics/chromosome_communities/Review1/paths_to_consider.txt'



#### Old

For each chromosome:
- top: annotation bars;
- middle: collapsed untangling output;
- bottom: average entropy.

```shell
Rscript /lizardfs/guarracino/chromosome_communities/scripts/figures/OldFig4_Annotation_CollapsedUntangle_AggregatedEntropy.R \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.support.dedup.eid0900.n1.nref1.tsv.gz \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.entropy.by_contig.eid0900.w50000.n1.nref1.tsv \
  /lizardfs/guarracino/chromosome_communities/data/annotation/ \
  ~/Figure4.pdf
  
Rscript /lizardfs/guarracino/chromosome_communities/scripts/figures/Fig4_Annotation_CollapsedUntangle_AggregatedEntropy.R \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.support.dedup.eid0900.n1.nref1.tsv.gz \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.entropy.by_contig.eid0900.w100000.n1.nref1.tsv \
  /lizardfs/guarracino/chromosome_communities/data/annotation/ \
  ~/Figure4.w100k.pdf
```

### Figure 5

For each chromosome:
- top: annotation bars;
- bottom: best 5 untangled hits for selected contigs.

```shell
Rscript /lizardfs/guarracino/chromosome_communities/scripts/figures/OldFig5_Annotation_Untangle5hits.R \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.ALL.e50000.m1000.grounded.pq_touching.reliable.tsv.gz \
  /lizardfs/guarracino/chromosome_communities/data/annotation/ \
  0.90 \
  ~/Figure5.pdf
```

## Supplementary figures

## Figure 6, 7, 8, 9, 10

For each chromosome:
- top: annotation bars;
- bottom: best untangled hit for all pq-contigs.

```shell
path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.ALL.e50000.m1000.grounded.pq_touching.reliable.tsv.gz

n=6
(seq 13 15; seq 21 22) | while read i; do
    echo "chr$i"
        
    Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation.R \
      $path_grounded_pq_touching_reliable_tsv_gz \
      0 25000000 \
      90 0.7 \
      0 \
      1 1 \
      $i \
      0.9 \
      <(zgrep '^HG002#1\|^HG002#2' -v $path_grounded_pq_touching_reliable_tsv_gz | sed '1d' | cut -f 1 | sort | uniq) \
      /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \
      ~/SupplementaryFigure${n}.pdf
      
    n=$((n+1))
done
```

## Figure 12 (and 13)

Average entropy across chrX and Y:

```shell
path_entropy_by_contig_tsv=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/entropy/chrSEX+refs.fa.gz.2ed2c67.04f1c29.22fc5c8.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.reliable.entropy.by_contig.eid0900.w50000.n1.nref1.tsv

Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_entropy_with_BED_annotation.R \
  $path_entropy_by_contig_tsv \
  0 155000000 \
  90 \
  'X' \
  /lizardfs/guarracino/chromosome_communities/data/chm13_hg002.PARs.bed \
  ~/chrSEX+refs.fa.gz.2ed2c67.04f1c29.22fc5c8.smooth.final.untangle.chm13#chrSEX.e50000.m1000.grounded.reliable.entropy.eid0900.w50000.n1.nref1.chrX.pdf
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_entropy_with_BED_annotation.R \
  $path_entropy_by_contig_tsv \
  0 63000000 \
  90 \
  'Y' \
  /lizardfs/guarracino/chromosome_communities/data/chm13_hg002.PARs.bed \
  ~/chrSEX+refs.fa.gz.2ed2c67.04f1c29.22fc5c8.smooth.final.untangle.chm13#chrSEX.e50000.m1000.grounded.reliable.entropy.eid0900.w50000.n1.nref1.chrY.pdf
```

## Figure 19 (20,21,22,23)

```shell
# --delta for white space between the pieces
pdfjam --delta '0 7' --no-landscape --nup 1x5 \
  chrSEX+refs.fa.gz.2ed2c67.04f1c29.22fc5c8.smooth.final.untangle.chm13#chrSEX.e50000.m1000.grounded.reliable.entropy.eid0900.w50000.n1.nref1.chrX.pdf \
  chrSEX+refs.fa.gz.2ed2c67.04f1c29.22fc5c8.smooth.final.untangle.chm13#chrSEX.e50000.m1000.grounded.reliable.entropy.eid0900.w50000.n1.nref1.chrY.pdf \
  --outfile output.pdf

pdfcrop --margins "1 1 1 1" output.pdf SupplementaryFigure12.chrSEX+refs.fa.gz.2ed2c67.04f1c29.22fc5c8.smooth.final.untangle.chm13#chrSEX.e50000.m1000.grounded.reliable.entropy.eid0900.w50000.n1.nref1.chrSEX.pdf
rm output.pdf
```

## Figure 14, 15, 16, 17, 18

For each chromosome:
- top: annotation bars;
- bottom: best 5 untangled hits for all pq-contigs.

```shell
path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.ALL.e50000.m1000.grounded.pq_touching.reliable.tsv.gz

n=14
(seq 13 15; seq 21 22) | while read i; do
    echo "chr$i"
    
    Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation.R \
      $path_grounded_pq_touching_reliable_tsv_gz \
      0 25000000 \
      90 0.8 \
      0.6 \
      5 1 \
      $i \
      0.9 \
      <(zgrep '^HG002#1\|^HG002#2' -v $path_grounded_pq_touching_reliable_tsv_gz | sed '1d' | cut -f 1 | sort | uniq) \
      /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \
      ~/SupplementaryFigure${n}.pdf
         
    n=$((n+1))
done
```


## Figure 19 (20,21,22,23)

```shell
# --delta for white space between the pieces
pdfjam --delta '0 7' --no-landscape --nup 1x5 \
  chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.entropy_match_order.eid0900.n5.chr13.pdf \
  chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.entropy_match_order.eid0900.n5.chr14.pdf \
  chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.entropy_match_order.eid0900.n5.chr15.pdf \
  chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.entropy_match_order.eid0900.n5.chr21.pdf \
  chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.entropy_match_order.eid0900.n5.chr22.pdf \
  --outfile output.pdf

pdfcrop --margins "1 1 1 1" output.pdf SupplementaryFigure19.chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.entropy_match_order.eid0900.n5.chrACRO.pdf
rm output.pdf
```

## Figure 24 (25)

```shell
# --delta for white space between the pieces
pdfjam --delta '0 7' --no-landscape --nup 1x5 \
  chrSEX+refs.fa.gz.2ed2c67.04f1c29.22fc5c8.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.entropy_match_order.eid0900.n2.chrX.pdf \
  chrSEX+refs.fa.gz.2ed2c67.04f1c29.22fc5c8.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.entropy_match_order.eid0900.n2.chrY.pdf \
  --outfile output.pdf

pdfcrop --margins "1 1 1 1" output.pdf SupplementaryFigure24.chrSEX+refs.fa.gz.2ed2c67.04f1c29.22fc5c8.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.entropy_match_order.eid0900.n2.chrSEX.pdf
rm output.pdf
```

## Figure 28, 29, 30, 31, 32

For each chromosome:
- top: annotation bars;
- middle: concordance by haplotype;
- bottom: best untangled hit for HG002 contigs.

```shell
Rscript /lizardfs/guarracino/chromosome_communities/scripts/figures/SuppFig_Annotation_Concordance_UntangleBestHit.R \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.concordance.by_haplotype.eid0900.n1.nref1.tsv \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.ALL.e50000.m1000.grounded.pq_touching.reliable.tsv.gz \
  0.9 \
  /lizardfs/guarracino/chromosome_communities/data/annotation/ \
  ~
```


## Figure XXX

For each chromosome, length distribution of the untangled query segments.

```shell
path_grounded_pq_touching_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.ALL.e50000.m1000.grounded.pq_touching.reliable.tsv.gz

e=50000

Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangled_segment_histogram.R \
  $path_grounded_pq_touching_reliable_ALL_tsv_gz \
  0 25000000 \
  60 15 \
  $(echo "$e + 15000" | bc) \
  1 1 \
  0.9 \
  <( zgrep '^chm13\|^grch38\|^HG002#1\|HG002#2\|^HG01978#MAT\|^HG01978#PAT\|bakeoff' $path_grounded_pq_touching_reliable_ALL_tsv_gz -v | sed '1d' | cut -f 1 | sort | uniq ) \
  ~/SuppFigure17.eid090.0_25Mbps.pdf
```


## Figure XXX

Untangling dot plot for all contigs:

```shell
f=/lizardfs/guarracino/chromosome_communities/untangle/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#ACRO.e50000.m1000.j0.n100.fixed.bed.gz
(zcat $f | head -n 1; \
zcat $f | awk '$8 == "-" { x=$6; $6=$5; $5=x; } { print }' | awk '$10 == 1') | tr ' ' '\t' | pigz -c > for_dot_plot.bed.gz

Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_dot_plots.R \
  for_dot_plot.bed.gz \
  0.9 \
  ~/SuppFigureXX.DotPlots.0_25Mbps.eid090.pdf
```