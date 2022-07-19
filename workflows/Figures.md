# Figures

## Main figures

### Figure 4

For each chromosome:
- top: annotation bars;
- middle: collapsed untangling output;
- bottom: average entropy.

```shell
Rscript /lizardfs/guarracino/chromosome_communities/scripts/figures/Fig4_Annotation_CollapsedUntangle_AggregatedEntropy.R \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.support.dedup.eid0900.n1.nref1.tsv.gz \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.entropy.by_contig.eid0900.n1.nref1.tsv \
  /lizardfs/guarracino/chromosome_communities/data/annotation/ \
  ~
```

### Figure 5

For each chromosome:
- top: annotation bars;
- bottom: best 5 untangled hits for selected contigs.

```shell 
Rscript /lizardfs/guarracino/chromosome_communities/scripts/figures/Fig5_Annotation_Untangle5hits.R \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.ALL.e50000.m1000.grounded.pq_touching.reliable.tsv.gz \
  /lizardfs/guarracino/chromosome_communities/data/annotation/ \
  0.90 \
  ~/Figure5.pdf
```

## Supplementary figures

## Figure 2, 3, 4, 5, 6

For each chromosome:
- top: annotation bars;
- bottom: best untangled hit for all pq-contigs.

```shell
path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.ALL.e50000.m1000.grounded.pq_touching.reliable.tsv.gz

n=2
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
      ~/SuppFigure${n}.pdf
      
    n=$((n+1))
done
```

## Figure 7, 8, 9, 10, 11

For each chromosome:
- top: annotation bars;
- bottom: best 5 untangled hits for all pq-contigs.

```shell
path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.ALL.e50000.m1000.grounded.pq_touching.reliable.tsv.gz

n=7
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
      ~/SuppFigure${n}.pdf
         
    n=$((n+1))
done
```

## Figure 12, 13, 14, 15, 16

For each chromosome:
- top: annotation bars;
- middle: concordance by haplotype;
- bottom: best untangled hit for HG002 contigs.

```shell
Rscript /lizardfs/guarracino/chromosome_communities/scripts/figures/SuppFig_Annotation_Concordance_UntangleBestHit.R \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.concordance.by_haplotype.n1.nref1.tsv \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.ALL.e50000.m1000.grounded.pq_touching.reliable.tsv.gz \
  /lizardfs/guarracino/chromosome_communities/data/annotation/ \
  ~
```

## Figure 17

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
  ~/SuppFigure17.pdf
```


## Figure XXX

Untangling dot plot for all contigs:

```shell
f=/lizardfs/guarracino/chromosome_communities/untangle/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#ACRO.e50000.m1000.j0.n100.fixed.bed.gz
(zcat $f | head -n 1; \
zcat $f | awk '$8 == "-" { x=$6; $6=$5; $5=x; } { print }' | awk '$10 == 1') | tr ' ' '\t' | pigz -c > for_dot_plot.bed.gz

Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_dot_plots.R \
  for_dot_plot.bed.gz \
  ~
```