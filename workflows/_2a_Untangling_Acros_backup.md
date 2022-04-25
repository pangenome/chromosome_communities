# Untangling

## Acrocentric chromosomes

### Collect contigs running from the p-arm to the q-arm of the acrocentric chromosomes

Prepare CHM13's acrocentric chromosomes:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/assemblies
cd /lizardfs/guarracino/chromosome_communities/assemblies

wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/CHM13/assemblies/chm13.draft_v1.1.fasta.gz
~/tools/fastix/target/release/fastix-331c1159ea16625ee79d1a82522e800c99206834 -p 'chm13#' <(zcat chm13.draft_v1.1.fasta.gz) |\
  bgzip -@ 48 -c >chm13.fa.gz
samtools faidx chm13.fa.gz
rm chm13.draft_v1.1.fasta.gz

(seq 13 15; seq 21 22) | while read f; do
  echo $chr$f
  samtools faidx chm13.fa.gz $(echo chm13#chr$f) | bgzip -@ 48 -c > chm13.chr$f.fa.gz && samtools faidx chm13.chr$f.fa.gz
done
```

Map acrocentric contigs against the acrocentric CHM13's chromosomes:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/mappings
cd /lizardfs/guarracino/chromosome_communities/mappings

for s in 300k 200k 150k 100k 50k 20k ; do
  for p in 98 95 90 85 80; do
    s_no_k=${s::-1}
    l_no_k=$(echo $s_no_k '*' 3 | bc)
    l=${l_no_k}k
    
    (seq 13 15; seq 21 22) | while read f; do
      if [[ ! -s /lizardfs/guarracino/chromosome_communities/mappings/chr$f.vs.chm13.chr$f.s$s.l$l.p$p.n1.paf ]]; then
        sbatch -p allnodes -c 24 --job-name AcroVsRef --wrap 'hostname; cd /scratch; \time -v /gnu/store/d6ajy0ajxkbb22gkgi32g5c4jy8rs8bv-wfmash-0.6.0+948f168-21/bin/wfmash /lizardfs/guarracino/chromosome_communities/assemblies/chm13.chr'$f'.fa.gz /lizardfs/erikg/HPRC/year1v2genbank/parts/chr'$f'.pan.fa -s '$s' -l '$l' -p '$p' -n 1 -t 24 -m > /lizardfs/guarracino/chromosome_communities/mappings/chr'$f'.vs.chm13.chr'$f'.s'$s'.l'$l'.p'$p'.n1.paf'
      fi;
    done
  done
done

# Debugging
#ls *pq_contigs.paf | while read f; do
#	~/git/wfmash/scripts/paf2dotplot png large /home/guarracino/Downloads/Pangenomics/chromosome_communities/mappings/$f
# 	mv out.png $f.png
#done
```


```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/pq_contigs
cd /lizardfs/guarracino/chromosome_communities/pq_contigs
```

Find the contigs which have mappings at least 1kbp-long in both the p-arm and the q-arm of the same chromosome:

```shell
for s in 300k 200k 150k 100k 50k 20k ; do
  for p in 98 95 90 85 80; do
    s_no_k=${s::-1}
    l_no_k=$(echo $s_no_k '*' 3 | bc)
    l=${l_no_k}k
    
    (seq 13 15; seq 21 22) | while read f; do      
      if [[ ! -s /lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.chr$f.s$s.l$l.p$p.n1.pq_contigs.txt ]]; then
        echo s$s l$l p$p chr$f
          
        # "p-touching contigs" intersected "q-touching contigs"
        comm -12 \
          <(bedtools intersect \
            -a <(cat /lizardfs/guarracino/chromosome_communities/mappings/chr$f.vs.chm13.chr$f.s$s.l$l.p$p.n1.paf | awk -v OFS='\t' '{print $6, $8, $9, $1, "", "+"}' | bedtools sort) \
            -b <(grep chr$f p_arms.bed | bedtools sort) | \
            awk '$3-$2+1>=1000' | \
            cut -f 4 | sort | uniq) \
          <(bedtools intersect \
            -a <(cat /lizardfs/guarracino/chromosome_communities/mappings/chr$f.vs.chm13.chr$f.s$s.l$l.p$p.n1.paf | awk -v OFS='\t' '{print $6, $8, $9, $1, "", "+"}' | bedtools sort) \
            -b <(grep chr$f q_arms.bed | bedtools sort) | \
            awk '$3-$2+1>=1000' | \
            cut -f 4 | sort | uniq) | \
            grep 'chm13#\|grch38#' -v > /lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.chr$f.s$s.l$l.p$p.n1.1kbps.pq_contigs.txt
      fi;
    done
  done
done
```

Take the union of pq-contigs found with all combinations of parameters:

```shell
(seq 13 15; seq 21 22) | while read f; do
  cat /lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.chr$f.s$s.l$l.p$p.n1.1kbps.pq_contigs.txt |\
   sort | uniq > /lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.1kbps.pq_contigs.union.txt
done

# Num. of pq-contigs
#(seq 13 15; seq 21 22) | while read f; do echo -n "$f -> "; cat /lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.1kbps.pq_contigs.union.txt | wc -l; done
```

Prepare the FASTA files with the pq-contigs:

```shell
cd /lizardfs/guarracino/chromosome_communities/pq_contigs

(seq 13 15; seq 21 22) | while read f; do
  path_pq_contigs_fa_gz=/lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.1kbps.pq_contigs.union.hg002prox.fa.gz
  if [[ ! -s $path_pq_contigs_fa_gz ]]; then
    echo chr$f

    path_pq_contigs_txt=/lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.1kbps.pq_contigs.union.txt
  
    cat \
      <(samtools faidx /lizardfs/erikg/HPRC/year1v2genbank/parts/chr$f.pan.fa $(cat $path_pq_contigs_txt)) \
      <(zcat /lizardfs/guarracino/chromosome_communities/assemblies/hg002.chr$f.prox.fa.gz) | \
      bgzip -@ 48 -c > $path_pq_contigs_fa_gz
    samtools faidx $path_pq_contigs_fa_gz
  fi;
done

# Num. of contigs
#(seq 13 15; seq 21 22) | while read f; do echo -n "$f -> "; cat /lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.1kbps.pq_contigs.union.hg002prox.fa.gz.fai | wc -l; done
```

Put all together with both the human references:

```shell
path_grch38_acro_fa_gz=/lizardfs/guarracino/chromosome_communities/pq_contigs/grch38.ACRO.fa.gz
samtools faidx /lizardfs/erikg/HPRC/year1v2genbank/assemblies/grch38.fa $(grep 'chr13\|chr14\|chr15\|chr21\|chr22' /lizardfs/erikg/HPRC/year1v2genbank/assemblies/grch38.fa.fai | grep '_' -v | cut -f 1) | \
  bgzip -@ 48 -c > $path_grch38_acro_fa_gz

zcat \
  /lizardfs/guarracino/chromosome_communities/assemblies/chm13.chr*.fa.gz \
  $path_grch38_acro_fa_gz \
  /lizardfs/guarracino/chromosome_communities/pq_contigs/chr*.vs.chm13.1kbps.pq_contigs.union.hg002prox.fa.gz | \
  
  bgzip -@ 48 -c > /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.fa.gz
samtools faidx /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.fa.gz
```

Put all together with both the human references plus all partitioned HG002 HiFi contigs:

```shell
zcat \
  /lizardfs/guarracino/chromosome_communities/assemblies/chm13.chr*.fa.gz \
  $path_grch38_acro_fa_gz \
  /lizardfs/guarracino/chromosome_communities/pq_contigs/chr*.vs.chm13.1kbps.pq_contigs.union.hg002prox.fa.gz \
  /lizardfs/guarracino/chromosome_communities/assemblies/hg002.chr*.hifi.fa.gz | \
  bgzip -@ 48 -c > /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.hg002hifi.fa.gz
samtools faidx /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.hg002hifi.fa.gz

# With `HG002#2#h2tg000134l         4689773  0  4689773  -  chm13#chr21           45090682   6659180   11348953  42408  4689773  10  id:f:90.4268`
zcat \
  /lizardfs/guarracino/chromosome_communities/assemblies/chm13.chr*.fa.gz \
  $path_grch38_acro_fa_gz \
  /lizardfs/guarracino/chromosome_communities/pq_contigs/chr*.vs.chm13.1kbps.pq_contigs.union.hg002prox.fa.gz \
  /lizardfs/guarracino/chromosome_communities/assemblies/hg002.chr*.hifi.fa.gz \
  <( samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/HG002.maternal.fa.gz HG002#2#h2tg000134l | bgzip -c ) | \
  bgzip -@ 48 -c > /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.hg002hifi.h2tg000134.fa.gz
samtools faidx /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.hg002hifi.h2tg000134.fa.gz
```

### Collect all acrocentric contigs

```shell
# Prepare sequence order, with all references on the top
grep chr /lizardfs/erikg/HPRC/year1v2genbank/parts/chrA.pan.fa.fai | cut -f 1 > sequence_order.txt
grep chr /lizardfs/erikg/HPRC/year1v2genbank/parts/chrA.pan.fa.fai -v | cut -f 1 >> sequence_order.txt

cat \
  <(zcat hg002.chr13.prox.fa.gz) \
  <(zcat hg002.chr14.prox.fa.gz) \
  <(zcat hg002.chr15.prox.fa.gz) \
  <(zcat hg002.chr21.prox.fa.gz) \
  <(zcat hg002.chr22.prox.fa.gz) \
  <(samtools faidx /lizardfs/erikg/HPRC/year1v2genbank/parts/chrA.pan.fa $(cat sequence_order.txt)) |\
   bgzip -@ 48 > chrA.pan+HG002chrAprox.fa.gz
samtools faidx chrA.pan+HG002chrAprox.fa.gz

rm sequence_order.txt
```

### Pangenome building

Apply `pggb` on the pq-contigs:
##################################################################### Add -F 0.01 to pggb
```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/graphs

num_of_haplotypes=$(cat /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.fa.gz.fai | wc -l)
sbatch -p highmem -c 48 --job-name acropggb --wrap 'hostname; cd /scratch && /gnu/store/swnkjnc9wj6i1cl9iqa79chnf40r1327-pggb-0.2.0+640bf6b-5/bin/pggb -i /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.fa.gz -o chrACRO+refs.1kbps.pq_contigs.union.hg002prox.s100k.l300k.p98.n'$num_of_haplotypes' -t 48 -s 100k -l 300k -p 98 -n '$num_of_haplotypes' -k 311 -G 13117,13219 -O 0.03 -T 48 -v -V chm13:#,grch38:#; mv /scratch/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.s100k.l300k.p98.n'$num_of_haplotypes' /lizardfs/guarracino/chromosome_communities/graphs';

num_of_haplotypes=$(cat /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.fa.gz.fai | wc -l)
num_of_haplotypes_plus_a_bit=$(echo "$num_of_haplotypes + 10" | bc) # 5 mat acros + 5 pat acros
sbatch -p highmem -c 48 --job-name acropggb --wrap 'hostname; cd /scratch && /gnu/store/swnkjnc9wj6i1cl9iqa79chnf40r1327-pggb-0.2.0+640bf6b-5/bin/pggb -i /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.hg002hifi.fa.gz -o chrACRO+refs.1kbps.pq_contigs.union.hg002prox.hg002hifi.s100k.l300k.p98.n'$num_of_haplotypes_plus_a_bit' -t 48 -s 100k -l 300k -p 98 -n '$num_of_haplotypes_plus_a_bit' -k 311 -G 13117,13219 -O 0.03 -T 48 -v -V chm13:#,grch38:#; mv /scratch/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.hg002hifi.s100k.l300k.p98.n'$num_of_haplotypes_plus_a_bit' /lizardfs/guarracino/chromosome_communities/graphs';

num_of_haplotypes=$(cat /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.fa.gz.fai | wc -l)
num_of_haplotypes_plus_a_bit=$(echo "$num_of_haplotypes + 10" | bc) # 5 mat acros + 5 pat acros
sbatch -p highmem -c 48 --job-name acropggb --wrap 'hostname; cd /scratch && /gnu/store/swnkjnc9wj6i1cl9iqa79chnf40r1327-pggb-0.2.0+640bf6b-5/bin/pggb -i /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.hg002hifi.h2tg000134.fa.gz -o chrACRO+refs.1kbps.pq_contigs.union.hg002prox.hg002hifi.h2tg000134.s100k.l300k.p98.n'$num_of_haplotypes_plus_a_bit' -t 48 -s 100k -l 300k -p 98 -n '$num_of_haplotypes_plus_a_bit' -k 311 -G 13117,13219 -O 0.03 -T 48 -v -V chm13:#,grch38:#; mv /scratch/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.hg002hifi.h2tg000134.s100k.l300k.p98.n'$num_of_haplotypes_plus_a_bit' /lizardfs/guarracino/chromosome_communities/graphs';


#num_of_haplotypes=$(cut -f 1,2 -d '#' /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.fai | sort | uniq | wc -l)
#sbatch -p highmem -c 48 --job-name acropggb --wrap 'hostname; cd /scratch && /gnu/store/swnkjnc9wj6i1cl9iqa79chnf40r1327-pggb-0.2.0+640bf6b-5/bin/pggb -i /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz -o chrACRO+refs.50kbps.pq_contigs.union.hg002prox.s100k.l300k.p98.n'$num_of_haplotypes' -t 48 -s 100k -l 300k -p 98 -n '$num_of_haplotypes' -k 311 -G 13117,13219 -O 0.03 -T 48 -v -V chm13:#,grch38:#; mv /scratch/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.s100k.l300k.p98.n'$num_of_haplotypes' /lizardfs/guarracino/chromosome_communities/graphs';
```

### Untangling

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/untangle

#path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.s100k.l300k.p98.n178/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.fa.gz.de31bcf.4030258.2385969.smooth.fix.og
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.hg002hifi.s100k.l300k.p98.n188/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.hg002hifi.fa.gz.4817bf7.4030258.57ed14c.smooth.fix.og

prefix=$(basename $path_input_og .og)
run_odgi=/home/guarracino/tools/odgi/bin/odgi-694948ccf31e7b565449cc056668e9dcc8cc0a3e

path_targets_txt=/lizardfs/guarracino/chromosome_communities/untangle/chm13.target_paths.txt
grep chm13 /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.fa.gz.fai | cut -f 1 > $path_targets_txt

# All references and emit cut points
for e in 50000 5000 100000; do
  for m in 1000 500 10000; do
    echo "-e $e -m $m"
    
    path_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.bed.gz
    path_cut_points_txt=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.cut_points.txt
    
    if [[ ! -s ${path_cut_points_txt} ]]; then
      sbatch -p workers -c 24 --job-name acrountangle --wrap '\time -v '$run_odgi' untangle -t 24 -P -i '$path_input_og' -R '$path_targets_txt' -e '$e' -m '$m' --cut-points-output '$path_cut_points_txt' -j 0 -n 1 | pigz -c > '$path_bed_gz';'
    fi;
  done
done

# Single reference by using the same cut points
for e in 50000 5000 100000; do
  for m in 1000 500 10000; do
    echo "-e $e -m $m"
    
    path_cut_points_txt=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.cut_points.txt
      
    cat $path_targets_txt | while read ref; do
      echo $ref
      
      path_ref_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$ref.e$e.m$m.j0.n100.bed.gz
      if [[ ! -s ${path_ref_bed_gz} ]]; then
        path_cut_points_txt=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.cut_points.txt
        
        sbatch -p workers -c 24 --job-name acrountangle --wrap '\time -v '$run_odgi' untangle -t 24 -P -i '$path_input_og' -r '$ref' -e '$e' -m '$m' --cut-points-input '$path_cut_points_txt' -j 0 -n 100 | pigz -c > '$path_ref_bed_gz';'
      fi;
    done
  done
done
```

Grounding (applying filters) and plotting:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/untangle/grounded

for e in 50000 ; do
  for m in 1000 ; do
    for j in 0 0.8 0.95; do
      j_str=$(echo $j | sed 's/\.//g')
      (seq 1 5; seq 10 10 50) | while read n; do 
        echo "-e $e -m $m -j $j -n $n"
          
        cat $path_targets_txt | while read ref; do
          echo -e "\t"$ref
            
          path_grounded_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.tsv.gz
            
          if [[ ! -s ${path_grounded_tsv_gz} ]]; then
            # Grounding
            ( echo query query.begin query.end target target.begin target.end jaccard strand self.coverage ref ref.begin ref.end | tr ' ' '\t'
              join \
                <(zcat /lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.bed.gz | awk '{ print $1"_"$2, $0 }' | tr ' ' '\t' | sort -k 1,1) \
                <(zcat /lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$ref.e$e.m$m.j0.n100.bed.gz | awk -v j=$j -v n=$n '{ if($7 >= j && $10 <= n) {print $1"_"$2, $0} }' | tr ' ' '\t' | sort -k 1,1) | \
              tr ' ' '\t' | grep -v '^#' | cut -f 2- | cut -f -9,14-16 ) | tr ' ' '\t' | pigz -c > ${path_grounded_tsv_gz}
          fi;
        done
      done
    done
  done
done

# Take pq-untangling contigs
for e in 50000 ; do
  for m in 1000 ; do
    for j in 0 0.8 0.95; do
      j_str=$(echo $j | sed 's/\.//g')
      (seq 1 5; seq 10 10 50) | while read n; do 
        echo "-e $e -m $m -j $j -n $n"
    
        cat $path_targets_txt | while read ref; do
          echo $ref
              
          path_grounded_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.tsv.gz
          path_grounded_pq_touching_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv
              
          if [[ ! -s ${path_grounded_pq_touching_tsv}.gz ]]; then
            # "p-touching contigs" intersected "q-touching contigs"
            comm -12 \
              <(bedtools intersect \
                -a <(zcat ${path_grounded_tsv_gz} | awk -v OFS="\t" '{print $10,$11,$12,$1, "", "+"}' | grep query -v | bedtools sort) \
                -b <(grep $ref /lizardfs/guarracino/chromosome_communities/pq_contigs/p_arms.bed | bedtools sort) | \
                #awk '$3-$2+1>=100000' | \
                cut -f 4 | sort | uniq) \
              <(bedtools intersect \
                -a <(zcat ${path_grounded_tsv_gz} | awk -v OFS="\t" '{print $10,$11,$12,$1, "", "+"}' | grep query -v | bedtools sort) \
                -b <(grep $ref /lizardfs/guarracino/chromosome_communities/pq_contigs/q_arms.bed | bedtools sort) | \
                #awk '$3-$2+1>=100000' | \
                cut -f 4 | sort | uniq) | \
                grep 'chm13#\|grch38#' -v > $ref.tmp.txt

            # Consider only chromosome-partitioned contigs
            comm -12 \
              <(sort $ref.tmp.txt) \
              <(cut -f 1 /lizardfs/guarracino/chromosome_communities/pq_contigs/`echo $ref | sed 's/chm13#//g'`.vs.chm13.1kbps.pq_contigs.union.hg002prox.fa.gz.fai | sort)  \
              > $ref.tmp2.txt
            rm $ref.tmp.txt
            
            ##########################################################################################
            # To plot the short (bad) HiFi-based HG002 contigs
            cat /lizardfs/guarracino/chromosome_communities/assemblies/hg002.`echo $ref | sed 's/chm13#//g'`.hifi.fa.gz.fai | awk '$2 > 300000' | cut -f 1 >> $ref.tmp2.txt       
            if [ $ref == "chm13#chr13" ]; then
#              echo "Remove wrongly partitioned contig: HG002#1#h1tg000013l"
#              grep 'HG002#1#h1tg000013l' -v $ref.tmp2.txt > $ref.tmp3.txt
#              echo "Add wrongly partitioned contig: HG002#1#h1tg000260l"
#              echo 'HG002#1#h1tg000260l' >> $ref.tmp3.txt

              echo "Remove wrongly partitioned contig: HG002#1#JAHKSE010000013.1"
              grep 'HG002#1#JAHKSE010000013.1' -v $ref.tmp2.txt > $ref.tmp3.txt
              echo "Add wrongly partitioned contig: HG002#1#JAHKSE010000214.1"
              echo 'HG002#1#JAHKSE010000214.1' >> $ref.tmp3.txt             
              
              rm $ref.tmp2.txt && mv $ref.tmp3.txt $ref.tmp2.txt
            fi
            if [ $ref == "chm13#chr21" ]; then
#              echo "Remove wrongly partitioned contig: HG002#1#h1tg000260l"
#              grep 'HG002#1#h1tg000260l' -v $ref.tmp2.txt > $ref.tmp3.txt
#              echo "Add wrongly partitioned contig: HG002#1#h1tg000013l"
#              echo 'HG002#1#h1tg000013l' >> $ref.tmp3.txt

              echo "Remove wrongly partitioned contig: HG002#1#JAHKSE010000214.1"
              grep 'HG002#1#JAHKSE010000214.1' -v $ref.tmp2.txt > $ref.tmp3.txt
              echo "Add wrongly partitioned contig: HG002#1#JAHKSE010000013.1"
              echo 'HG002#1#JAHKSE010000013.1' >> $ref.tmp3.txt 
              
              rm $ref.tmp2.txt && mv $ref.tmp3.txt $ref.tmp2.txt
            fi              
            ##########################################################################################
          
            #### Put header (with a new 'target' column), take intersection, and re-add other acrocentric references
            cat \
              <(zcat ${path_grounded_tsv_gz} | head -n 1 | awk -v OFS='\t' '{print $0, "grounded.target"}') \
              <(zgrep ${path_grounded_tsv_gz} -f $ref.tmp2.txt | awk -v OFS='\t' -v ref=$ref '{print $0, ref}') \
              <(zgrep '^grch38\|^chm13' ${path_grounded_tsv_gz}  | awk -v OFS='\t' -v ref=$ref '{print $0, ref}') \
                > ${path_grounded_pq_touching_tsv}
            rm $ref.tmp2.txt
                  
            # Add annotation
            zgrep rDNA /lizardfs/guarracino/chromosome_communities/data/chm13.CentromeresAndTelomeres.CenSatAnnotation.txt.gz | \
              sed 's/chr/chm13#chr/g' | \
              grep $ref | \
              awk -v OFS='\t' -v ref=$ref '{print $1"#"$4,".",".",$1,".",".",".",".","1","ref",$2,$3, ref}' >> ${path_grounded_pq_touching_tsv}
            cat /lizardfs/guarracino/chromosome_communities/data/chm13.centromeres.approximate.bed | \
              bedtools merge | \
              grep 'chr13\|chr14\|chr15\|chr21\|chr22' | \
              sed 's/chr/chm13#chr/g'  | \
              grep $ref | \
              awk -v OFS='\t' -v ref=$ref '{print $1"#centromere",".",".",$1,".",".",".",".","1","ref",$2,$3, ref}' >> ${path_grounded_pq_touching_tsv}
                  
            pigz ${path_grounded_pq_touching_tsv}
          fi;
        done
      done
    done
  done
done
```

Remove unreliable regions:

```shell
# Remove unreliable regions:
for e in 50000 ; do
  for m in 1000 ; do
    for j in 0.8 0.95; do
      j_str=$(echo $j | sed 's/\.//g')
      (seq 1 5; seq 10 10 50) | while read n; do 
        echo "-e $e -m $m -j $j -n $n"
    
        cat $path_targets_txt | while read ref; do
          echo $ref

          path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz
          path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.reliable.tsv.gz

          if [[ ! -s ${path_grounded_pq_touching_reliable_tsv_gz} ]]; then
            # Skip verkko's HG002 contigs because we don't have unreliable regions for it
            zcat $path_grounded_pq_touching_tsv_gz | sed '1d' | grep 'HG002#MAT\|HG002#PAT' -v | cut -f 1 | sort | uniq | while read CONTIG; do
              SAMPLE=$( echo $CONTIG | cut -f 1 -d '#')
            
              path_unreliable_bed=/lizardfs/erikg/HPRC/year1v2genbank/annotations/unreliable/$SAMPLE.hifi.flagger_final.bed
              if [[ -s $path_unreliable_bed ]]; then
                #echo $CONTIG "--->" $SAMPLE
                
                zgrep "^$CONTIG" $path_grounded_pq_touching_tsv_gz | sort -k 2n | awk -v OFS='\t' '{print($1,$2,$3,$4"_"$5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11"_"$12"_"$13)}' > x.bed
                grep $CONTIG $path_unreliable_bed > y.bed
                # -A: remove entire feature if any overlap
                bedtools subtract -a x.bed -b y.bed -A |\
                  awk -v OFS='\t' '{split($4, a, "_"); print($1,$2,$3,a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10])}' >> x.tsv
              else
                zgrep "^$CONTIG" $path_grounded_pq_touching_tsv_gz | sort -k 2n >> x.tsv
              fi
            done
            
            # Re-take verkko's HG002 untangled regions
            zcat $path_grounded_pq_touching_tsv_gz | sed '1d' | grep 'HG002#MAT\|HG002#PAT' >> x.tsv
            
            cat \
              <( zcat $path_grounded_pq_touching_tsv_gz | head -n 1 ) \
               x.tsv | pigz -c > $path_grounded_pq_touching_reliable_tsv_gz
            rm x.tsv
          fi;
        done
      done
    done
  done
done

```

Merged plots:

```shell
# Dependencies
#guix install r
#guix install r-ggplot2
#guix install r-ggforce
#guix install poppler # For pdfunite

mkdir -p /lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/

for e in 50000 ; do
  for m in 1000 ; do
    for j in 0.8 0.95; do
        j_str=$(echo $j | sed 's/\.//g')
        (seq 1 5; seq 10 10 50) | while read n; do 
            echo "-e $e -m $m -j $j -n $n"
    
            path_grounded_pq_touching_reliable_all_chromosomes_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chrACRO.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.reliable.tsv.gz
            if [[ ! -s ${path_grounded_pq_touching_reliable_all_chromosomes_tsv_gz} ]]; then
                # Merge single reference results
                cat \
                  <(zcat /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.*.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.reliable.tsv.gz | head -n 1) \
                  <(zcat /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.*.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.reliable.tsv.gz | grep query -v) |\
                  pigz -c > x.tsv.gz
                # Rename after to avoid getting itself with the previous '*' expansion
                mv x.tsv.gz ${path_grounded_pq_touching_reliable_all_chromosomes_tsv_gz}
            fi;

            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_all.R ${path_grounded_pq_touching_reliable_all_chromosomes_tsv_gz} "-e $e -m $m -j $j -n $n" 0 25000000 120 200
        done
      done
    done
done

mv /lizardfs/guarracino/chromosome_communities/untangle/grounded/*.pdf /lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/

# Merge
for e in 50000 ; do
  for m in 1000 ; do
    for j in 0.8 0.95; do
        echo "-e $e -m $m -j $j"
        
        j_str=$(echo $j | sed 's/\.//g')
        PDFs=$((seq 1 5; seq 10 10 50) | while read n; do \
          echo /lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/$prefix.untangle.chm13#chrACRO.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.reliable.tsv.gz.pdf
        done | tr '\n' ' ')
        #echo $PDFs
        
        # Too slow
        #gs -sDEVICE=pdfwrite \
        #  -dNOPAUSE -dBATCH -dSAFER \
        #  -sOutputFile=/lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/$prefix.untangle.chm13#chrACRO.e$e.m$m.j${j_str}.grounded.pq_touching.tsv.gz.pdf \
        #  $PDFs
          
        /gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite $PDFs /lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/$prefix.untangle.chm13#chrACRO.e$e.m$m.j${j_str}.merged.grounded.pq_touching.reliable.tsv.gz.pdf
      done
    done
done
```

Single plots (not used):

```shell
#for e in 5000 50000 100000; do
#  for m in 500 1000 10000; do
#      for j in 0 0.8; do
#        j_str=$(echo $j | sed 's/\.//g')
#        (seq 1 5; seq 10 10 50) | while read n; do 
#            echo "-e $e -m $m -j $j -n $n"
#    
#            grep chm13 /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.fa.gz.fai | cut -f 1 | while read ref; do
#              echo $ref
#              
#              path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz
#            
#              Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle.R ${path_grounded_pq_touching_tsv_gz} "$ref -e $e -m $m -j $j -n $n"
#            done
#        done
#      done
#    done
#done
```


Add Mobin's annotations:
Take only reliable blocks (flagged with "Hh" or "Hc", https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/README.md#components)

```shell
e=50000
m=1000
j=0.8
j_str=$(echo $j | sed 's/\.//g')
n=1
ref=chm13#chr13

# Brutal: remove untangled regions if they overlap with unreliable regions
touch xyz.tsv
cat $path_targets_txt | while read ref; do
    echo $ref
    path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz
    
    touch z.tsv
    zcat $path_grounded_pq_touching_tsv_gz | sed '1d' | grep 'HG002#MAT\|HG002#PAT' -v | cut -f 1 | sort | uniq | while read CONTIG; do
      SAMPLE=$( echo $CONTIG | cut -f 1 -d '#')
    
      path_unreliable_bed=/lizardfs/erikg/HPRC/year1v2genbank/annotations/unreliable/$SAMPLE.hifi.flagger_final.bed
      if [[ -s $path_unreliable_bed ]]; then
        echo $CONTIG "--->" $SAMPLE
        
        zgrep "^$CONTIG" $path_grounded_pq_touching_tsv_gz | sort -k 2n | awk -v OFS='\t' '{print($1,$2,$3,$4"_"$5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11"_"$12"_"$13)}' > x.bed
        grep $CONTIG $path_unreliable_bed > y.bed
        # -A: remove entire feature if any overlap
        bedtools subtract -a x.bed -b y.bed -A |\
          awk -v OFS='\t' '{split($4, a, "_"); print($1,$2,$3,a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10])}' >> xyz.tsv
      else
        zgrep "^$CONTIG" $path_grounded_pq_touching_tsv_gz | sort -k 2n >> xyz.tsv
      fi
    done
    
    zcat $path_grounded_pq_touching_tsv_gz | sed '1d' | grep 'HG002#MAT\|HG002#PAT' >> xyz.tsv
done

path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz
cat \
  <( zcat $path_grounded_pq_touching_tsv_gz | head -n 1 ) \
   xyz.tsv | pigz -c > xyz.tsv.gz
rm xyz.tsv

Rscript scripts/plot_untangle_all.R xyz.tsv.gz "-e 50000 -m 1000 -j 0.8 -n 1" 0 25000000 360 200


# Plot additional tracks for the annotations
touch xyz.tsv
cat $path_targets_txt | while read ref; do
    echo $ref
    path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz
    
    touch z.tsv
    zcat $path_grounded_pq_touching_tsv_gz | sed '1d' | grep 'HG002#MAT\|HG002#PAT' -v | cut -f 1 | sort | uniq | while read CONTIG; do
      SAMPLE=$( echo $CONTIG | cut -f 1 -d '#')
    
      path_unreliable_bed=/lizardfs/erikg/HPRC/year1v2genbank/annotations/unreliable/$SAMPLE.hifi.flagger_final.bed
      if [[ -s $path_unreliable_bed ]]; then
        echo $CONTIG "--->" $SAMPLE
        
        zgrep "^$CONTIG" $path_grounded_pq_touching_tsv_gz | sort -k 2n | awk -v OFS='\t' '{print($1,$2,$3,$4"_"$5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11"_"$12"_"$13)}' > x.bed
        grep $CONTIG $path_unreliable_bed > y.bed
        
        # wao: write the original A and B entries plus the number of base pairs of overlap between the two features.
        bedtools intersect -a x.bed -b y.bed -wo >> z.tsv
      fi
    done
    
    cat \
      <( zcat $path_grounded_pq_touching_tsv_gz | sed '1d' ) \
      <( python3 scripts/get_annotation_track.py z.tsv ) | tr ' ' '\t' >> xyz.tsv

    rm z.tsv
done

path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz
cat \
  <( zcat $path_grounded_pq_touching_tsv_gz | head -n 1 ) \
   xyz.tsv | pigz -c > xyz.annot.tsv.gz
rm xyz.tsv

Rscript scripts/plot_untangle_all.R xyz.annot.tsv.gz "-e 50000 -m 1000 -j 0.8 -n 1" 0 25000000 120 200
```

## Variant calling

Call variants in a haploid setting:

```shell
# pq-contigs
path_input_gfa=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.s100k.l300k.p98.n84/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.gfa
path_chm13_vcf_gz=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.s100k.l300k.p98.n84/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.vcf.gz
#path_grch38_vcf_gz=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.s100k.l300k.p98.n84/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.grch38.vcf.gz
sbatch -p workers -c 48 --job-name vgchm13 --wrap '\time -v vg deconstruct -P chm13 -H '?' -e -a -t 48 '$path_input_gfa' | bgzip -@ 48 -c > '$path_chm13_vcf_gz' && tabix '$path_chm13_vcf_gz
#sbatch -p workers -c 48 --job-name vggrch38 --wrap '\time -v vg deconstruct -P grch38 -H '?' -e -a -t 48 '$path_input_gfa' | bgzip -@ 48 -c > '$path_grch38_vcf_gz' && tabix '$path_grch38_vcf_gz

# acro-contigs
path_input_gfa=/lizardfs/guarracino/chromosome_communities/graphs/chrA.pan+HG002chrAprox.s100k.l300k.p98.n100/chrA.pan+HG002chrAprox.fa.gz.952ce3b.4030258.0e65f78.smooth.fix.gfa
path_input_sed_gfa=/lizardfs/guarracino/chromosome_communities/graphs/chrA.pan+HG002chrAprox.s100k.l300k.p98.n100/chrA.pan+HG002chrAprox.fa.gz.952ce3b.4030258.0e65f78.smooth.fix.sed.gfa
path_chm13_vcf_gz=/lizardfs/guarracino/chromosome_communities/graphs/chrA.pan+HG002chrAprox.s100k.l300k.p98.n100/chrA.pan+HG002chrAprox.fa.gz.952ce3b.4030258.0e65f78.smooth.fix.sed.c1000.chm13.vcf.gz
#path_grch38_vcf_gz=/lizardfs/guarracino/chromosome_communities/graphs/chrA.pan+HG002chrAprox.s100k.l300k.p98.n100/chrA.pan+HG002chrAprox.fa.gz.952ce3b.4030258.0e65f78.smooth.fix.sed.grch38.vcf.gz

# sed replaces only the first instance on a line by default (without the /g modifier)
# To have names like NA21309-1#1#JAHEPC010000450.1 and call haploid genotypes with -H
sed 's/#/-/' $path_input_gfa | sed 's/#/#1#/' > $path_input_sed_gfa

sbatch -p headnode -c 48 --job-name vgchm13 --wrap '\time -v vg deconstruct -P chm13 -H "#" -e -a -c 1000 -t 48 '$path_input_sed_gfa' | bgzip -@ 48 -c > '$path_chm13_vcf_gz' && tabix '$path_chm13_vcf_gz
#sbatch -p highmem -c 48 --job-name vggrch38 --wrap '\time -v vg deconstruct -P grch38 -H "#" -e -c 1000 -a -t 48 '$path_input_sed_gfa' | bgzip -@ 48 -c > '$path_grch38_vcf_gz' && tabix '$path_grch38_vcf_gz
```