# Figures

## Main figures

### Figure 4

For each chromosome:
- top: annotation bars;
- middle: collapsed untangling output;
- bottom: average entropy.

```shell
Rscript /lizardfs/guarracino/chromosome_communities/scripts/figures/Fig4_Annotation_CollapsedUntangle_AggregatedEntropy.R \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.support.dedup.n1.nref1.tsv.gz \
  /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.untangle.chm13#chrACRO.e50000.m1000.grounded.pq_touching.reliable.entropy.by_contig.n1.nref1.tsv \
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
  ~
done
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
      91 0.8 \
      0 \
      1 1 \
      $i \
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
      91 0.8 \
      0.8 \
      5 1 \
      $i \
      <(zgrep '^HG002#1\|^HG002#2' -v $path_grounded_pq_touching_reliable_tsv_gz | sed '1d' | cut -f 1 | sort | uniq) \
      /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \
      ~/SuppFigure${n}.pdf
         
    n=$((n+1))
done
```
