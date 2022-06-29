# Figures

## Figure 4

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

## Figure 5

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
