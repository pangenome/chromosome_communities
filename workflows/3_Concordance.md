# Concordance

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance
cd /lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance
```

Compute the concordance between verkko's HG002 and HiFi's HG002:

```shell
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.og
prefix=$(basename $path_input_og .og)

mkdir -p /lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/

for e in 50000 ; do
  for m in 1000 ; do
    echo "-e $e -m $m"
  
    # By haplotype
    path_concordance_by_haplotype_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.concordance.by_haplotype.tsv
    if [[ ! -s ${path_concordance_by_haplotype_tsv} ]]; then
      python3 /lizardfs/guarracino/chromosome_communities/scripts/concordance.by_haplotype.py $path_grounded_pq_touching_reliable_all_chromosomes_tsv_gz chm13#ACRO.len.tsv 1 1 > $path_concordance_by_haplotype_tsv
    fi
    
    PREFIX=$(basename $path_concordance_by_haplotype_tsv .tsv);
    (seq 13 15; seq 21 22) | while read i; do
      Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_concordance_tile_with_annotation.by_haplotype.R \
        $path_concordance_by_haplotype_tsv \
        0 25000000 \
        100 1 \
        $i \
        /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \
        /lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/$PREFIX.chr${i}.pdf
    done
    
    # Merge chromosomes's PDF files
    /gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite \
      /lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/$PREFIX.chr*.pdf \
      /lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/$PREFIX.merged.pdf
  

    # By contig
    path_grounded_pq_touching_reliable_all_chromosomes_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.tsv.gz

    path_concordance_by_contig_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.concordance.by_contig.tsv
    if [[ ! -s ${path_concordance_by_contig_tsv} ]]; then
      python3 /lizardfs/guarracino/chromosome_communities/scripts/concordance.by_contig.py $path_grounded_pq_touching_reliable_all_chromosomes_tsv_gz chm13#ACRO.len.tsv 1 1 > $path_concordance_by_contig_tsv
    fi

    PREFIX=$(basename $path_concordance_by_contig_tsv .tsv);
    (seq 13 15; seq 21 22) | while read i; do
      Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_concordance_tile_with_annotation.by_contig.R \
        $path_concordance_by_contig_tsv \
        0 25000000 \
        80 1.2 \
        $i \
        /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \
        /lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/$PREFIX.chr${i}.pdf
    done
    
    # Merge chromosomes's PDF files
    /gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite \
      /lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/$PREFIX.chr*.pdf \
      /lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/$PREFIX.merged.pdf
  done
done
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
