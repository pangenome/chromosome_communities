### Entropy

Compute the entropy for each contig and for each chromosome (by averaging each window across contigs):

```shell
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.og
prefix=$(basename $path_input_og .og)

mkdir -p /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/

#guix install r-dplyr
for e in 50000; do
  for m in 1000 2000; do
    for refn in 1 10; do
      echo "-e $e -m $m -refn $refn"

      path_entropy_by_contig_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.entropy.by_contig.n1.nref${refn}.tsv
      if [[ ! -s ${path_entropy_by_contig_tsv} ]]; then
        path_grounded_pq_touching_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.ALL.e$e.m$m.grounded.pq_touching.reliable.tsv.gz
        python3 /lizardfs/guarracino/chromosome_communities/scripts/entropy.by_contig.py \
          $path_grounded_pq_touching_reliable_ALL_tsv_gz \
          chm13#ACRO.len.tsv \
          1 $refn > $path_entropy_by_contig_tsv
      fi

      PREFIX=$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.entropy.n1.nref${refn}
      (seq 13 15; seq 21 22) | while read i; do
        Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_entropy_tile_with_annotation.by_contig.R \
          $path_entropy_by_contig_tsv \
          0 25000000 \
          98 1.2 \
          $i \
          /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \
          /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$PREFIX.chr${i}
      done
    
      # Merge chromosomes's PDF files
      /gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite \
        /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$PREFIX.chr*.by_contig.pdf \
        /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$PREFIX.by_contig.merged.pdf
      rm /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$PREFIX.chr*.by_contig.pdf
      
      /gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite \
        /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$PREFIX.chr*.by_chromosome.pdf \
        /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$PREFIX.by_chromosome.merged.pdf
      rm /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$PREFIX.chr*.by_chromosome.pdf
    done
  done
done
```


### Concordance

Compute the concordance between verkko's HG002 and HiFi's HG002:

```shell
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.og
prefix=$(basename $path_input_og .og)

mkdir -p /lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/

# By haplotype
for e in 50000; do
  for m in 1000 2000; do
    for refn in 1 10; do
      echo "-e $e -m $m -refn $refn"
    
      path_concordance_by_haplotype_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.concordance.by_haplotype.n1.nref${refn}.tsv
      if [[ ! -s ${path_concordance_by_haplotype_tsv} ]]; then
        path_grounded_pq_touching_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.ALL.e$e.m$m.grounded.pq_touching.reliable.tsv.gz
        python3 /lizardfs/guarracino/chromosome_communities/scripts/concordance.by_haplotype.py \
          $path_grounded_pq_touching_reliable_ALL_tsv_gz \
          chm13#ACRO.len.tsv \
          1 $refn > $path_concordance_by_haplotype_tsv
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
      rm /lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/$PREFIX.chr*.pdf
    done
  done
done

# By contig
for e in 50000; do
  for m in 1000 2000; do
    for refn in 1 10; do
      echo "-e $e -m $m -refn $refn"

      path_concordance_by_contig_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.concordance.by_contig.n1.nref${refn}.tsv
      if [[ ! -s ${path_concordance_by_contig_tsv} ]]; then
        path_grounded_pq_touching_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.ALL.e$e.m$m.grounded.pq_touching.reliable.tsv.gz
        python3 /lizardfs/guarracino/chromosome_communities/scripts/concordance.by_contig.py \
          $path_grounded_pq_touching_reliable_ALL_tsv_gz \
          chm13#ACRO.len.tsv \
          1 $refn > $path_concordance_by_contig_tsv
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
      rm /lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/$PREFIX.chr*.pdf
    done
  done
done
```