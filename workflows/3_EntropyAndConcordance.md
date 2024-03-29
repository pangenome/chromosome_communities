### Entropy

Compute the entropy for each contig and for each chromosome (by averaging each window across contigs) 
by considering HiFi-only contigs anchored to the q-arms (so no HG002-HiFi-only) and HG002-verkko:


#### Acrocentric chromosomes

```shell
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.og
prefix=$(basename $path_input_og .og)

mkdir -p /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/

#guix install r-dplyr
for e in 50000; do
  for m in 1000; do
    for eid in 0.900 0.950 0.975 0.995 1.000; do
      eid_str=$(echo $eid | sed 's/\.//g')
      
      for w in 50000 100000; do
        for refn in 1; do   
          path_entropy_by_contig_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.entropy.by_contig.eid${eid_str}.w${w}.n1.nref${refn}.tsv
          if [[ ! -s ${path_entropy_by_contig_tsv} ]]; then
            echo "-e $e -m $m $eid -w $w -refn $refn"
            path_grounded_pq_touching_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.ALL.e$e.m$m.grounded.pq_touching.reliable.tsv.gz
            python3 /lizardfs/guarracino/chromosome_communities/scripts/entropy.by_contig.py \
              $path_grounded_pq_touching_reliable_ALL_tsv_gz \
              /lizardfs/guarracino/chromosome_communities/chm13#ACRO.len.tsv \
              $w \
              1 $refn \
              $eid \
              <( zgrep '^chm13\|^grch38\|^HG002#1\|HG002#2\|^HG01978#MAT\|^HG01978#PAT\|bakeoff' $path_grounded_pq_touching_reliable_ALL_tsv_gz -v | sed '1d' | cut -f 1 | sort | uniq ) \
              > $path_entropy_by_contig_tsv
          fi
        
          # Plot
          PREFIX=$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.entropy.eid${eid_str}.w${w}.n1.nref${refn}
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
  done
done
```


#### Sex chromosomes

```shell
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrSEX+refs.s50k.l250k.p98.n102/chrSEX+refs.fa.gz.2ed2c67.04f1c29.22fc5c8.smooth.final.og
prefix=$(basename $path_input_og .og)

mkdir -p /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/entropy/

#guix install r-dplyr
for e in 50000; do
  for m in 1000; do
    for eid in 0.900; do
      eid_str=$(echo $eid | sed 's/\.//g')
      
      for w in 50000 100000; do
        for refn in 1; do
          echo "-e $e -m $m $eid -refn $refn"
    
          path_entropy_by_contig_tsv=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/entropy/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.reliable.entropy.by_contig.eid${eid_str}.w${w}.n1.nref${refn}.tsv
          if [[ ! -s ${path_entropy_by_contig_tsv} ]]; then
            path_grounded_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.ALL.e$e.m$m.grounded.reliable.tsv.gz
            python3 /lizardfs/guarracino/chromosome_communities/scripts/entropy.by_contig.py \
              $path_grounded_reliable_ALL_tsv_gz \
              /lizardfs/guarracino/chromosome_communities/chm13#SEX.len.tsv \
              $w \
              1 $refn \
              $eid \
              <( zgrep '^chm13\|^grch38\|^HG002#1\|HG002#2\|^HG01978#MAT\|^HG01978#PAT\|bakeoff' $path_grounded_reliable_ALL_tsv_gz -v | sed '1d' | cut -f 1 | sort | uniq ) \
              > $path_entropy_by_contig_tsv
          fi
        
          # Plot
          PREFIX=$prefix.untangle.chm13#chrSEX.e$e.m$m.grounded.reliable.entropy.eid${eid_str}.w${w}.n1.nref${refn}
          (echo X; echo Y) | while read i; do
            xmax=155000000
            if [[ $i == "Y" ]]; then
              xmax=63000000
            fi
            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_entropy_with_BED_annotation.R \
              $path_entropy_by_contig_tsv \
              0 $xmax \
              90 \
              $i \
              /lizardfs/guarracino/chromosome_communities/data/chm13_hg002.PARs.bed \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/entropy/$PREFIX.chr${i}.pdf
          done
          
          # Merge chromosomes's PDF files
          #/gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite \
          #  /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/entropy/$PREFIX.chr*.pdf \
          #  /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/entropy/$PREFIX.merged.pdf
          #rm /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/entropy/$PREFIX.chr*.pdf
        done
      done
    done
  done
done
```


### Entropy match order

Compute the entropy for each reference position for all samples untangling in that position
by considering HiFi-only contigs anchored to the q-arms (so no HG002-HiFi-only) and HG002-verkko.

#### Acrocentric chromosomes

```shell
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.og
prefix=$(basename $path_input_og .og)


n=5

for e in 50000; do
  for m in 1000; do
    path_grounded_pq_touching_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.ALL.e$e.m$m.grounded.pq_touching.reliable.tsv.gz
    PREFIX=$(basename $path_grounded_pq_touching_reliable_ALL_tsv_gz .tsv.gz)
      
    zgrep '^chm13\|^grch38\|^HG002#1\|HG002#2\|^HG01978#MAT\|^HG01978#PAT\|bakeoff' $path_grounded_pq_touching_reliable_ALL_tsv_gz -v | sed '1d' | cut -f 1 | sort | uniq \
      > /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$PREFIX.query_to_consider.txt
    
    for eid in 0.900 0.950 0.975 0.995 1.000; do
      eid_str=$(echo $eid | sed 's/\.//g')
      
      echo "-e $e -m $m $eid"

      path_entropy_match_order_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.entropy_match_order.eid${eid_str}.n${n}.tsv
      if [[ ! -s ${path_entropy_match_order_tsv} ]]; then
        python3 /lizardfs/guarracino/chromosome_communities/scripts/entropy_match_orders.py \
          $path_grounded_pq_touching_reliable_ALL_tsv_gz \
          /lizardfs/guarracino/chromosome_communities/chm13#ACRO.len.tsv \
          $n $eid \
          /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$PREFIX.query_to_consider.txt \
          > $path_entropy_match_order_tsv
      fi
    done
  done
done
```

Compute the Pseudo-Homologous Regions (PHRs):

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/PHRs

e=50000
m=1000
eid=0.900
n=5
eid_str=$(echo $eid | sed 's/\.//g')
path_entropy_match_order_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.entropy_match_order.eid${eid_str}.n${n}.tsv
      
cat $path_entropy_match_order_tsv | sed '1d' | awk '$4 > 0 && $5 > 0' | \
  bedtools merge -i - -d 30000 -c 4,5 -o mean | \
  awk '$3 - $2 > 30000' | sed 's/chm13#//' > /lizardfs/guarracino/chromosome_communities/PHRs/chrACRO_7-Dec-22_PHRs.bed # Supplementary Table

# To check the sizes
f=/lizardfs/guarracino/chromosome_communities/PHRs/chrACRO_7-Dec-22_PHRs.bed
cut -f 1 $f | sort | uniq | while read x; do
  echo $f $x
  grep $x $f | awk '{l+=$3-$2;}END{print(l/1000000)}'
done
```

Plots:

```shell
#guix install r-dplyr

n=5
for e in 50000; do
  for m in 1000; do
    for eid in 0.900; do
      eid_str=$(echo $eid | sed 's/\.//g')
      
      path_entropy_match_order_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.entropy_match_order.eid${eid_str}.n${n}.tsv
      PREFIX=$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.entropy_match_order.eid${eid_str}.n${n}
      
      echo "-e $e -m $m $eid"

      (seq 13 15; seq 21 22) | while read i; do              
        Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_entropy_match_order_with_annotation.R \
          $path_entropy_match_order_tsv \
          0 25000000 \
          90 \
          $i \
          /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \
          /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$PREFIX.chr${i}.pdf
      done
      
      # Merge chromosomes's PDF files
      PREFIX=$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.entropy_match_order.eid${eid_str}.n${n}
      #/gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite \
      #  /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$PREFIX.chr*.pdf \
      #  /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$PREFIX.merged.pdf
      #rm /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$PREFIX.chr*.pdf
    
      (seq 13 15; seq 21 22) | while read i; do    
        Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_entropy_match_order_with_BED_annotation.R \
          $path_entropy_match_order_tsv \
          0 25000000 \
          90 \
          $i \
          /lizardfs/guarracino/chromosome_communities/PHRs/chrACRO_7-Dec-22_PHRs.bed \
          /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$PREFIX.chr${i}.PHRs.pdf
      done
        
      # Merge chromosomes's PDF files
      #PREFIX=$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.entropy_match_order.eid${eid_str}.n${n}
      #/gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite \
      #  /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$PREFIX.chr*.PHRs.pdf \
      #  /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$PREFIX.PHRs.merged.pdf
      #rm /lizardfs/guarracino/chromosome_communities/untangle/grounded/entropy/$PREFIX.chr*.PHRs.pdf
     done
  done
done
```

#### Sex chromosomes

```shell
cd /lizardfs/guarracino/chromosome_communities/

mkdir -p /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/entropy/

path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrSEX+refs.s50k.l250k.p98.n102/chrSEX+refs.fa.gz.2ed2c67.04f1c29.22fc5c8.smooth.final.og
prefix=$(basename $path_input_og .og)


n=2

for e in 50000; do
  for m in 1000; do
    path_grounded_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.ALL.e$e.m$m.grounded.reliable.tsv.gz
    PREFIX=$(basename $path_grounded_pq_touching_reliable_ALL_tsv_gz .tsv.gz)
      
    zgrep '^chm13\|^grch38\|^HG002#1\|HG002#2\|^HG01978#MAT\|^HG01978#PAT\|bakeoff' $path_grounded_reliable_ALL_tsv_gz -v | sed '1d' | cut -f 1 | sort | uniq \
      > /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/entropy/$PREFIX.query_to_consider.txt
    
    for eid in 0.900; do
      eid_str=$(echo $eid | sed 's/\.//g')
      
      echo "-e $e -m $m $eid"

      path_entropy_match_order_tsv=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/entropy/$prefix.untangle.chm13#chrSEX.e$e.m$m.grounded.pq_touching.reliable.entropy_match_order.eid${eid_str}.n${n}.tsv
      if [[ ! -s ${path_entropy_match_order_tsv} ]]; then
        python3 /lizardfs/guarracino/chromosome_communities/scripts/entropy_match_orders.py \
          $path_grounded_reliable_ALL_tsv_gz \
          /lizardfs/guarracino/chromosome_communities/chm13#SEX.len.tsv \
          $n $eid \
          /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/entropy/$PREFIX.query_to_consider.txt \
          > $path_entropy_match_order_tsv
      fi
    done
  done
done
```

Compute the Pseudo-Homologous Regions (PHRs):

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/PHRs

e=50000
m=1000
eid=0.900
n=2
eid_str=$(echo $eid | sed 's/\.//g')
path_entropy_match_order_tsv=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/entropy/$prefix.untangle.chm13#chrSEX.e$e.m$m.grounded.pq_touching.reliable.entropy_match_order.eid${eid_str}.n${n}.tsv

cat $path_entropy_match_order_tsv | sed '1d' | awk '$4 > 0 && $5 > 0' | \
  bedtools merge -i - -d 30000 -c 4,5 -o mean | \
  awk '$3 - $2 > 30000' | sed 's/chm13#//' > /lizardfs/guarracino/chromosome_communities/PHRs/chrSEX_7-Dec-22_PHRs.bed # Supplementary Table

# To check the sizes
f=/lizardfs/guarracino/chromosome_communities/PHRs/chrSEX_7-Dec-22_PHRs.bed 
cut -f 1 $f | sort | uniq | while read x; do
  echo $f $x
  grep $x $f | awk '{l+=$3-$2;}END{print(l/1000000)}'
done
```

Plots:

```shell
n=2

#guix install r-dplyr
for e in 50000; do
  for m in 1000; do
    for eid in 0.900; do
      eid_str=$(echo $eid | sed 's/\.//g')

      path_entropy_match_order_tsv=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/entropy/$prefix.untangle.chm13#chrSEX.e$e.m$m.grounded.pq_touching.reliable.entropy_match_order.eid${eid_str}.n${n}.tsv
      PREFIX=$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.entropy_match_order.eid${eid_str}.n${n}
        
      echo "-e $e -m $m $eid"

      (echo X; echo Y) | while read i; do
        xmax=155000000
        if [[ $i == "Y" ]]; then
          xmax=63000000
        fi
        
        Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_entropy_match_order_with_BED_annotation.R \
          $path_entropy_match_order_tsv \
          0 $xmax \
          90 \
          $i \
          /lizardfs/guarracino/chromosome_communities/data/chm13_hg002.PARs.bed \
          /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/entropy/$PREFIX.chr${i}.pdf
      done
      
      # Merge chromosomes's PDF files
      #/gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite \
      #  /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/entropy/$PREFIX.chr*.pdf \
      #  /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/entropy/$PREFIX.merged.pdf
      #rm /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/entropy/$PREFIX.chr*.pdf
    done
  done
done
# For Supplementari Figures 24 and 25
```


### Concordance

Compute the concordance between verkko's HG002 and HiFi's HG002:

```shell
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.og
prefix=$(basename $path_input_og .og)

mkdir -p /lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/

# By haplotype
for e in 50000; do
  for m in 1000; do
    for eid in 0.900 0.950 0.975 0.995 1.000; do
      eid_str=$(echo $eid | sed 's/\.//g')
      
      for refn in 1; do
        echo "-e $e -m $m $eid -refn $refn"
        
        path_concordance_by_haplotype_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.concordance.by_haplotype.eid${eid_str}.n1.nref${refn}.tsv
        PREFIX=$(basename $path_concordance_by_haplotype_tsv .tsv);
        if [[ ! -s ${path_concordance_by_haplotype_tsv} ]]; then
          path_grounded_pq_touching_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.ALL.e$e.m$m.grounded.pq_touching.reliable.tsv.gz
          python3 /lizardfs/guarracino/chromosome_communities/scripts/concordance.by_haplotype.py \
            $path_grounded_pq_touching_reliable_ALL_tsv_gz \
            /lizardfs/guarracino/chromosome_communities/chm13#ACRO.len.tsv \
            1 $refn $eid > $path_concordance_by_haplotype_tsv
        fi
            
        
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
        
        
        # Statistics whole chromosome
        python3 /lizardfs/guarracino/chromosome_communities/scripts/get_concordance_ratio.py $path_concordance_by_haplotype_tsv | \
          awk -v OFS='\t' '{print( (NR==1) ? "region" : "whole",$0)}' > /lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/$PREFIX.stats.tsv
        
        # Statistics p arm
        python3 /lizardfs/guarracino/chromosome_communities/scripts/get_concordance_ratio.py <(echo 'fake_header'; bedtools intersect \
          -a <( cat $path_concordance_by_haplotype_tsv | sed '1d' | awk -v OFS='\t' '{print($1,$2,$3,$4"_"$5)}' ) \
          -b <( sed 's/^chr/chm13#chr/g' /lizardfs/guarracino/chromosome_communities/data/annotation/chm13.p_arms.approximate.acros.bed) |\
           awk -v OFS='\t' '{split($4,a,/_/); print($1,$2,$3,a[1],a[2])}') | sed '1d' | \
           awk -v OFS='\t' '{print("p.arm",$0)}' >> /lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/$PREFIX.stats.tsv
        
        # Statistics q arm
        python3 /lizardfs/guarracino/chromosome_communities/scripts/get_concordance_ratio.py <(echo 'fake_header'; bedtools intersect \
          -a <( cat $path_concordance_by_haplotype_tsv | sed '1d' | awk -v OFS='\t' '{print($1,$2,$3,$4"_"$5)}' ) \
          -b <( sed 's/^chr/chm13#chr/g' /lizardfs/guarracino/chromosome_communities/data/annotation/chm13.q_arms.approximate.acros.bed) |\
           awk -v OFS='\t' '{split($4,a,/_/); print($1,$2,$3,a[1],a[2])}') | sed '1d' | \
           awk -v OFS='\t' '{print("q.arm",$0)}' >> /lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/$PREFIX.stats.tsv
      done
    done
  done
done

#> require(tidyverse)
#> read.delim('xxx.stats.tsv.tsv') %>% filter(region=="whole") %>% summarize(1-sum(num.bp.not.ok)/(sum(num.bp.ok)+sum(num.bp.not.ok)))
#1 0.9942085
#> read.delim('xxx.stats.tsv.tsv') %>% filter(region=="q.arm") %>% summarize(1-sum(num.bp.not.ok)/(sum(num.bp.ok)+sum(num.bp.not.ok)))
#1 0.9993142
#> read.delim('xxx.stats.tsv.tsv') %>% filter(region=="p.arm") %>% summarize(1-sum(num.bp.not.ok)/(sum(num.bp.ok)+sum(num.bp.not.ok)))
#1 0.8745065


# By contig
for e in 50000; do
  for m in 1000; do
    for eid in 0.900 0.950 0.975 0.995 1.000; do
      eid_str=$(echo $eid | sed 's/\.//g')
      
      for refn in 1; do
        echo "-e $e -m $m $eid -refn $refn"

        path_concordance_by_contig_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/concordance/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.concordance.by_contig.eid${eid_str}.n1.nref${refn}.tsv
        if [[ ! -s ${path_concordance_by_contig_tsv} ]]; then
          path_grounded_pq_touching_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.ALL.e$e.m$m.grounded.pq_touching.reliable.tsv.gz
          python3 /lizardfs/guarracino/chromosome_communities/scripts/concordance.by_contig.py \
            $path_grounded_pq_touching_reliable_ALL_tsv_gz \
            /lizardfs/guarracino/chromosome_communities/chm13#ACRO.len.tsv \
            1 $refn $eid > $path_concordance_by_contig_tsv
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
done
```
