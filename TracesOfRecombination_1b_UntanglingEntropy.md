# Traces of recombination

Take acrocentric chromosome lengths:

```shell
grep '^chm13' assemblies/chrA.pan+HG002chrAprox.fa.gz.fai | cut -f 1,2 > chm13#ACRO.len.tsv

python3 scripts/entropy.py scripts/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.fa.gz.de31bcf.4030258.2385969.smooth.fix.untangle.chm13#chrACRO.e50000.m10000.j08.n1.grounded.pq_touching.tsv.gz scripts/chm13#ACRO.len.tsv > x.tsv

Rscript scripts/plot_entropy.R x.tsv 'e50000.m10000.j08.n1 - Window 50kbp'
```