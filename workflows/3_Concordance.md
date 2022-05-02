# Concordance

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance
cd /lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance
```

Take acrocentric chromosome lengths:

```shell
grep '^chm13' /lizardfs/guarracino/chromosome_communities/assemblies/chrA.pan+HG002chrAprox.fa.gz.fai | cut -f 1,2 > chm13#ACRO.len.tsv
```

Compute the concordance between verkko's HG002 and HiFi's HG002:

```shell
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.og
prefix=$(basename $path_input_og .og)

for e in 50000 ; do
  for m in 1000 ; do
    for j in 0.8 0.95; do
      j_str=$(echo $j | sed 's/\.//g')
      n=1
      echo "-e $e -m $m -j $j -n $n"
      
      path_grounded_pq_touching_reliable_all_chromosomes_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chrACRO.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.reliable.tsv.gz

      path_concordance_by_contig_tsv=$prefix.untangle.chm13#chrACRO.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.reliable.concordance.by_contig.tsv
      if [[ ! -s ${path_concordance_by_contig_tsv} ]]; then
        python3 /lizardfs/guarracino/chromosome_communities/scripts/concordance.by_contig.py $path_grounded_pq_touching_reliable_all_chromosomes_tsv_gz chm13#ACRO.len.tsv > $path_concordance_by_contig_tsv
      fi
      Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_concordance_tile.by_contig.R $path_concordance_by_contig_tsv "-e $e -m $m -j $j -n $n" 25000000 70

      #path_concordance_tsv=$prefix.untangle.chm13#chrACRO.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.reliable.concordance.tsv
      #if [[ ! -s ${path_concordance_tsv} ]]; then
      #  python3 /lizardfs/guarracino/chromosome_communities/scripts/concordance.by_haplotype.py $path_grounded_pq_touching_reliable_all_chromosomes_tsv_gz chm13#ACRO.len.tsv > $path_concordance_tsv
      #fi
      #Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_concordance_tile.R $prefix.untangle.chm13#chrACRO.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.reliable.concordance.tsv "-e $e -m $m -j $j -n $n" 25000000 25
    done
  done
done


python3 scripts/concordance.by_contig.py chrACRO+refs.1kbps.pq_contigs.union.hg002prox.hg002hifi.fa.gz.4817bf7.4030258.57ed14c.smooth.fix.untangle.chm13#chrACRO.e50000.m1000.j08.n1.grounded.pq_touching.reliable.tsv.gz chm13#ACRO.len.tsv > concordance.by_contig.tsv

Rscript scripts/plot_concordance_tile.by_contig.R concordance.by_contig.tsv 'e50000.m1000.j08.n1' 25000000 100
```

```shell
zgrep '^HG002#MAT' chrACRO+refs.1kbps.pq_contigs.union.hg002prox.hg002hifi.fa.gz.4817bf7.4030258.57ed14c.smooth.fix.untangle.chm13#chrACRO.e50000.m1000.j08.n1.grounded.pq_touching.reliable.tsv.gz |\
 awk -v OFS='\t' '{print($13,$11,$12,$4)}' > HG002.maternal.untangle.bed

zgrep '^HG002#2#JAHKSD010000014.1' chrACRO+refs.1kbps.pq_contigs.union.hg002prox.hg002hifi.fa.gz.4817bf7.4030258.57ed14c.smooth.fix.untangle.chm13#chrACRO.e50000.m1000.j08.n1.grounded.pq_touching.reliable.tsv.gz |\
 awk -v OFS='\t' '{print($13,$11,$12,$4)}' > HG002#2#JAHKSD010000014.1.untangle.bed
```

Compute the entropy for each query:

```shell
for e in 50000 ; do
  for m in 1000 ; do
    for j in 0 0.8 0.95; do
      j_str=$(echo $j | sed 's/\.//g')
      n=1
      
      path_grounded_pq_touching_all_chromosomes_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chrACRO.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz
      python3 scripts/entropy.py $path_grounded_pq_touching_all_chromosomes_tsv_gz chm13#ACRO.len.tsv > entropy.tsv
        
      Rscript scripts/plot_entropy_tile.R entropy.tsv 'e50000.m1000.j08.n1 - Window 50kbp' 25000000 200
    done
  done
done
```
