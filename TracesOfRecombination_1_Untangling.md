# Traces of recombination

## Tools

```shell
mkdir -p ~/tools $$ cd ~/tools

git clone --recursive https://github.com/ekg/fastix.git
cd fastix
git checkout 331c1159ea16625ee79d1a82522e800c99206834
cargo build --release
mv target/release/fastix target/release/fastix-331c1159ea16625ee79d1a82522e800c99206834
cd ..

git clone --recursive https://github.com/pangenome/odgi.git
cd odgi
git checkout 358537678174adc975415a488b458725d8a213be
cmake -H. -Bbuild && cmake --build build -- -j 48
mv bin/odgi bin/odgi-358537678174adc975415a488b458725d8a213be
cd ..

(echo pggb wfmash seqwish smoothxg odgi gfaffix | tr ' ' '\n') | while read tool; do ls -l $(which $tool); done | cut -f 13 -d ' '
/gnu/store/swnkjnc9wj6i1cl9iqa79chnf40r1327-pggb-0.2.0+640bf6b-5/bin/pggb
/gnu/store/d6ajy0ajxkbb22gkgi32g5c4jy8rs8bv-wfmash-0.6.0+948f168-21/bin/wfmash
/gnu/store/v1q6ja2fy3c7fcz7bnr6k3x54amynhyp-seqwish-0.7.3+51ee550-1/bin/seqwish
/gnu/store/p8wflvwwxiih4z9ds5hfkin8mjld6qw2-smoothxg-0.6.0+0f383e5-10/bin/smoothxg
/gnu/store/2ln6zv8mk6aqqzcib39rgi11x2hn7mv9-odgi-0.6.2+9e9c481-13/bin/odgi
/gnu/store/ccw48k7h8v1brz3ap0sj3bcwvvmk6xra-gfaffix-0.1.2.2/bin/gfaffix

vg
vg: variation graph tool, version v1.36.0 "Cibottola"
```

## Preparation

Clone the repository:

```shell
cd /lizardfs/guarracino/
git clone --recursive https://github.com/pangenome/chromosome_communities.git
```

## Get coverage analysis BED files

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/coverage_analysis_y1_genbank
cd /lizardfs/guarracino/chromosome_communities/coverage_analysis_y1_genbank

prefix=https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/e9ad8022-1b30-11ec-ab04-0a13c5208311--COVERAGE_ANALYSIS_Y1_GENBANK/FLAGGER/JAN_09_2022/FINAL_HIFI_BASED/FLAGGER_HIFI_ASM_BEDS/

cut -f 1 /lizardfs/erikg/HPRC/year1v2genbank/parts/HPRCy1.pan.fa.gz.fai | cut -f 1 -d '#' | grep 'chm13\|grch38' -v | sort | uniq | while read sample; do
  wget -c $prefix$sample.hifi.flagger_final.bed
done
```

## Collect contigs running from the p-arm to the q-arm of the acrocentric chromosomes

Prepare CHM13's acrocentric chromosomes:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/assemblies
cd /lizardfs/guarracino/chromosome_communities/assemblies

wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/CHM13/assemblies/chm13.draft_v1.1.fasta.gz
~/tools/fastix/target/release/fastix-331c1159ea16625ee79d1a82522e800c99206834 -p 'chm13#' <(zcat chm13.draft_v1.1.fasta.gz) | bgzip -@ 48 -c >chm13.fa.gz
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

Prepare the BED files for p-arms and q-arms (-/+ 1Mbps):

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/pq_contigs
cd /lizardfs/guarracino/chromosome_communities/pq_contigs

bedtools slop \
    -i <(grep acen /lizardfs/guarracino/chromosome_communities/data/chm13.CytoBandIdeo.v2.txt | bedtools merge | sed 's/chr/chm13#chr/g' | sort) \
    -g <(cut -f 1,2 /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa.fai | sort) \
    -b 1000000 | \
  sort | \
  bedtools complement \
    -i - \
    -g <(cut -f 1,2 /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa.fai | sort) | \
  grep 'chr13\|chr14\|chr15\|chr21\|chr22' > tmp.bed
  
# Take odd rows
sed -n 1~2p tmp.bed > /lizardfs/guarracino/chromosome_communities/pq_contigs/p_arms.bed
# Take even rows
sed -n 2~2p tmp.bed > /lizardfs/guarracino/chromosome_communities/pq_contigs/q_arms.bed

rm tmp.bed
```

Find the contigs which have mappings at least 100kbp-long in both the p-arm and the q-arm of the same chromosome:

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
              awk '$3-$2+1>=100000' | \
              cut -f 4 | sort | uniq) \
            <(bedtools intersect \
              -a <(cat /lizardfs/guarracino/chromosome_communities/mappings/chr$f.vs.chm13.chr$f.s$s.l$l.p$p.n1.paf | awk -v OFS='\t' '{print $6, $8, $9, $1, "", "+"}' | bedtools sort) \
              -b <(grep chr$f q_arms.bed | bedtools sort) | \
              awk '$3-$2+1>=100000' | \
              cut -f 4 | sort | uniq) | \
              grep 'chm13#\|grch38#' -v > /lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.chr$f.s$s.l$l.p$p.n1.100kbps.pq_contigs.txt
      fi;
    done
  done
done
```

Take the union of pq-contigs found with all combinations of parameters:

```shell
(seq 13 15; seq 21 22) | while read f; do
  cat /lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.chr$f.s$s.l$l.p$p.n1.100kbps.pq_contigs.txt | sort | uniq > /lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.100kbps.pq_contigs.union.txt
done
```

Prepare the FASTA files with the pq-contigs:

```shell
(seq 13 15; seq 21 22) | while read f; do
  path_pq_contigs_fa_gz=/lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.100kbps.pq_contigs.union.fa.gz
  if [[ ! -s $path_pq_contigs_fa_gz ]]; then
    echo chr$f

    path_pq_contigs_txt=/lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.100kbps.pq_contigs.union.txt
  
    samtools faidx /lizardfs/erikg/HPRC/year1v2genbank/parts/chr$f.pan.fa $(cat $path_pq_contigs_txt) | \
      bgzip -@ 48 -c > $path_pq_contigs_fa_gz
    samtools faidx $path_pq_contigs_fa_gz
  fi;
done
```

Put all together with the references:

```shell
path_grch38_acro_fa_gz=/lizardfs/guarracino/chromosome_communities/pq_contigs/grch38.ACRO.fa.gz
samtools faidx /lizardfs/erikg/HPRC/year1v2genbank/assemblies/grch38.fa $(grep 'chr13\|chr14\|chr15\|chr21\|chr22' /lizardfs/erikg/HPRC/year1v2genbank/assemblies/grch38.fa.fai | grep '_' -v | cut -f 1) | \
  bgzip -@ 48 -c > $path_grch38_acro_fa_gz

zcat \
  /lizardfs/guarracino/chromosome_communities/assemblies/chm13.chr*.fa.gz \
  $path_grch38_acro_fa_gz \
  /lizardfs/guarracino/chromosome_communities/pq_contigs/chr*.vs.chm13.100kbps.pq_contigs.union.fa.gz | \
  bgzip -@ 48 -c > /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.100kbps.pq_contigs.union.fa.gz
samtools faidx /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.100kbps.pq_contigs.union.fa.gz
```

## Pangenome building

Apply `pggb` on the pq-contigs:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/graphs

num_of_haplotypes=$(cut -f 1,2 -d '#' /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.fai | sort | uniq | wc -l)
sbatch -p highmem -c 48 --job-name acropggb --wrap 'hostname; cd /scratch && /gnu/store/swnkjnc9wj6i1cl9iqa79chnf40r1327-pggb-0.2.0+640bf6b-5/bin/pggb -i /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.100kbps.pq_contigs.union.fa.gz -o chrACRO+refs.100kbps.pq_contigs.union.s100k.l300k.p98.n'$num_of_haplotypes' -t 48 -s 100k -l 300k -p 98 -n '$num_of_haplotypes' -k 311 -G 13117,13219 -O 0.03 -T 48 -v -V chm13:#,grch38:#; mv /scratch/chrACRO+refs.100kbps.pq_contigs.union.s100k.l300k.p98.n'$num_of_haplotypes' /lizardfs/guarracino/chromosome_communities/graphs';
```

## Untangling

```shell
# WIP
e=5000
m=1000
j=0.8
n=1

path_chr6_smooth=chr6.pan.fa.a2fb268.4030258.6a1ecc2.smooth.og
path_fasta=${path_chr6_smooth}.fasta
path_untangle_all_paf_gz=chr6.pan.fa.a2fb268.4030258.6a1ecc2.smooth.mhc.both.e$e.m$m.paf.gz
path_untangle_single_paf_gz=chr6.pan.fa.a2fb268.4030258.6a1ecc2.smooth.mhc.chm13.e$e.m$m.n$n.j$j.paf.gz
path_untangle_grounded_tsv_gz=chr6.mhc.untangle.chm13.e$e.m$m.n$n.j$j.grounded.tsv.gz
path_untangle_grounded_paf_gz=chr6.mhc.untangle.chm13.e$e.m$m.n$n.j$j.grounded.paf.gz
path_untangle_grounded_wfmash_paf_gz=chr6.mhc.untangle.chm13.e$e.m$m.n$n.j$j.grounded.wfmash.paf.gz

odgi paths -i ${path_chr6_smooth} -f > ${path_fasta}
samtools faidx ${path_fasta}

odgi untangle -i ${path_chr6_smooth} -t 16 -P -R <(echo chm13#chr6:28874656-33821293; echo grch38#chr6:28999849-34000009) -e $e -m $m -d cut_points.txt -p | pigz -c > ${path_untangle_all_paf_gz}
odgi untangle -i ${path_chr6_smooth} -t 16 -P -r chm13#chr6:28874656-33821293 -n $n -j $j -e $e -m $m -c cut_points.txt -p | pigz -c > ${path_untangle_single_paf_gz}

( echo query query.length query.begin query.end strand target target.length target.begin target.end num.res alignment.length map.q id jaccard self.coverage nb ref ref.length ref.begin ref.end | tr ' ' '\t' ; \
  join \
    <(zcat ${path_untangle_all_paf_gz} | awk '{ print $1"_"$3, $0 }' | tr ' ' '\t' | sort  -k 1,1) \
    <(zcat ${path_untangle_single_paf_gz} | awk -v OFS="\t" -v j=$j -v n=$n '{gsub("jc:f:", "", $14); gsub("nb:i:", "", $16); if($14 >= j && $16 <= n) {print $1"_"$3, $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,"jc:f:"$14,$15,"nb:i:"$16}}' | tr ' ' '\t' | sort  -k 1,1) | \
  tr ' ' '\t' | cut -f 2- | cut -f -16,23,24,25,26) | 
  tr ' ' '\t' | pigz -c > ${path_untangle_grounded_tsv_gz}

zgrep query -v ${path_untangle_grounded_tsv_gz} | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$17,$18,$19,$20,$10,$11,$12,$13,$14,$15,$16}' | pigz -c > ${path_untangle_grounded_paf_gz}
#zcat chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.untangle.chm13#chr13.e50000.m1000.n1.j08.grounded.paf.gz | -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,"id:f:99.99999",$14,$15,$16}' > chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.untangle.chm13#chr13.e50000.m1000.n1.j08.grounded.9999999.paf.gz
wfmash ${path_fasta} ${path_fasta} -i ${path_untangle_grounded_paf_gz} -t 6 | pigz -c > ${path_untangle_grounded_wfmash_paf_gz}

~/git/rustybam/target/release/rustybam stats -p <(zcat ${path_untangle_grounded_wfmash_paf_gz}) > chr6.mhc.untangle.chm13.e$e.m$m.n$n.j$j.grounded.wfmash.stats.bed
```

```
# WIP
# Saturate the alignments (-n 200)
sbatch -p highmem -c 48 --job-name acropggb --wrap 'hostname; cd /scratch && /gnu/store/swnkjnc9wj6i1cl9iqa79chnf40r1327-pggb-0.2.0+640bf6b-5/bin/pggb -i /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.100kbps.pq_contigs.union.fa.gz -o chrACRO+refs.100kbps.pq_contigs.union.s100k.l300k.p98.n200 -t 48 -s 100k -l 300k -p 98 -n 200 -k 311 -G 13117,13219 -O 0.03 -T 48 -v -V chm13:#,grch38:#; mv /scratch/chrACRO+refs.100kbps.pq_contigs.union.s100k.l300k.p98.n200 /lizardfs/guarracino/chromosome_communities/graphs';

############################
#!/bin/bash
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrA.pan.s100k.l300k.p97.n200/chrA.pan.fa.4825454.4030258.0517e7e.smooth.fix.og
prefix=$(basename $path_input_og .og)

for e in 50000; do
  for m in 1000; do
    echo "-e $e -m $m"
    
    path_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.bed.gz
    if [[ ! -s ${path_paf_gz} ]]; then
        path_cut_points_txt=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.cut_points.txt
        
        \time -v ~/tools/odgi/bin/odgi-358537678174adc975415a488b458725d8a213be untangle -t 48 -P \
          -i $path_input_og \
          -R <(grep chm13 /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.fai | cut -f 1) \
          -e $e -m $m \
          --cut-points-output ${path_cut_points_txt} | \
          pigz -c \
          > ${path_bed_gz}
    fi;
  done
done
############################

for ref in chm13#chr13 chm13#chr14 chm13#chr15 chm13#chr21 chm13#chr22; do
    echo $ref
    
    sbatch -p headnode -c 5 --job-name acrotangle untangle.sh $ref
done

untangle.sh############################
#!/bin/bash

ref=$1
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrA.pan.s100k.l300k.p97.n200/chrA.pan.fa.4825454.4030258.0517e7e.smooth.fix.og
prefix=$(basename $path_input_og .og)

for e in 50000; do
  for m in 1000; do
    echo "-e $e -m $m"
    
    path_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.bedpaf.gz
    if [[ ! -s ${path_paf_gz} ]]; then
        path_cut_points_txt=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.cut_points.txt

        echo $ref
        \time -v ~/tools/odgi/bin/odgi-358537678174adc975415a488b458725d8a213be untangle -t 5 -P \
          -i $path_input_og \
          -r $ref \
          -e $e -m $m \
          --cut-points-input ${path_cut_points_txt} \
          -j 0 -n 100 | pigz -c > /lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$ref.e$e.m$m.j0.n100.bed.gz
    fi;
  done
done
############################


# Grounding (applying filters)
for e in 50000; do
  for m in 1000; do
      for j in 0 0.8; do
        j_str=$(echo $j | sed 's/\.//g')
        (seq 10 10; seq 15 5 100) | while read n; do 
            echo "-e $e -m $m -j $j -n $n"
            
            for ref in chm13#chr13 chm13#chr14 chm13#chr15 chm13#chr21 chm13#chr22; do
              echo -e "\t"$ref
              
              path_grounded_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.n$n.j${j_str}.grounded.tsv.gz
              
              if [[ ! -s ${path_grounded_tsv_gz} ]]; then
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
for e in 50000; do
  for m in 1000; do
      for j in 0.8; do
        j_str=$(echo $j | sed 's/\.//g')
        (seq 10 10) | while read n; do 
            echo "-e $e -m $m -j $j -n $n"
    
            for ref in chm13#chr13 chm13#chr14 chm13#chr15 chm13#chr21 chm13#chr22; do
              echo $ref
              
              path_grounded_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.n$n.j${j_str}.grounded.tsv.gz
              path_grounded_pq_touching_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv
              
              if [[ ! -s ${path_grounded_pq_touching_tsv}.gz ]]; then
                  # "p-touching contigs" intersected "q-touching contigs"
                  comm -12 \
                    <(bedtools intersect \
                      -a <(zcat ${path_grounded_tsv_gz} | awk -v OFS="\t" '{print $10,$11,$12,$1, "", "+"}' | grep query -v | bedtools sort) \
                      -b <(cat /lizardfs/guarracino/chromosome_communities/pq_contigs/p_arms.bed | bedtools sort) | \
                      #awk '$3-$2+1>=100000' | \
                      cut -f 4 | sort | uniq) \
                    <(bedtools intersect \
                      -a <(zcat ${path_grounded_tsv_gz} | awk -v OFS="\t" '{print $10,$11,$12,$1, "", "+"}' | grep query -v | bedtools sort) \
                      -b <(cat /lizardfs/guarracino/chromosome_communities/pq_contigs/q_arms.bed | bedtools sort) | \
                      #awk '$3-$2+1>=100000' | \
                      cut -f 4 | sort | uniq) | \
                      grep 'chm13#\|grch38#' -v > $ref.tmp.txt
                
                  # Put header (with a new 'target' column), take intersection, and re-add other acrocentric references
                  cat \
                    <(zcat ${path_grounded_tsv_gz} | head -n 1 | awk -v OFS='\t' '{print $0, "grounded.target"}') \
                    <(zgrep ${path_grounded_tsv_gz} -f $ref.tmp.txt | awk -v OFS='\t' -v ref=$ref '{print $0, ref}') \
                    <(zgrep '^grch38\|^chm13' ${path_grounded_tsv_gz}  | awk -v OFS='\t' -v ref=$ref '{print $0, ref}') \
                     > ${path_grounded_pq_touching_tsv}
                
                  # Add annotation
                  zgrep rDNA /lizardfs/guarracino/chromosome_communities/data/chm13.CentromeresAndTelomeres.CenSatAnnotation.txt.gz | \
                    sed 's/chr/chm13#chr/g' | \
                    grep $ref | \
                    awk -v OFS='\t' -v ref=$ref '{print $1"#"$4,".",".",$1,".",".",".",".",".",".",$2,$3, ref}' >> ${path_grounded_pq_touching_tsv}
                  grep acen /lizardfs/guarracino/chromosome_communities/data/chm13.CytoBandIdeo.v2.txt | \
                    bedtools merge | \
                    grep 'chr13\|chr14\|chr15\|chr21\|chr22' | \
                    sed 's/chr/chm13#chr/g'  | \
                    grep $ref | \
                    awk -v OFS='\t' -v ref=$ref '{print $1"#centromere",".",".",$1,".",".",".",".",".",".",$2,$3, ref}' >> ${path_grounded_pq_touching_tsv}
                  
                  pigz ${path_grounded_pq_touching_tsv}
              fi;
            done
        done
      done
    done
done


# Merged plots
for e in 50000; do
  for m in 1000; do
      for j in 0 0.8; do
        j_str=$(echo $j | sed 's/\.//g')
        (seq 10 10) | while read n; do 
            echo "-e $e -m $m -j $j -n $n"
    
            path_grounded_pq_touching_all_chromosomes_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chrACRO.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz
            if [[ ! -s ${path_grounded_pq_touching_all_chromosomes_tsv_gz} ]]; then
                # Merge all acros together
                cat \
                  <(zcat /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.*.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz | head -n 1) \
                  <(zcat /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.*.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz | grep query -v) |\
                  pigz -c > x.tsv.gz
                mv x.tsv.gz ${path_grounded_pq_touching_all_chromosomes_tsv_gz}
            fi;

            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_all.R ${path_grounded_pq_touching_all_chromosomes_tsv_gz} "-e $e -m $m -j $j -n $n"
        done
      done
    done
done

mv /lizardfs/guarracino/chromosome_communities/untangle/grounded/*.pdf /lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/

`````

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/untangle

#!/bin/bash
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.100kbps.pq_contigs.union.s100k.l300k.p98.n97/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.og
prefix=$(basename $path_input_og .og)

for e in 5000 10000 20000 30000 40000 50000 100000; do
  for m in 1000 10000; do
    echo "-e $e -m $m"
    
    path_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.bedpaf.gz
    if [[ ! -s ${path_paf_gz} ]]; then
        path_cut_points_txt=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.cut_points.txt
        
        \time -v ~/tools/odgi/bin/odgi-358537678174adc975415a488b458725d8a213be untangle -t 48 -P \
          -i $path_input_og \
          -R <(grep chm13 /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.fai | cut -f 1) \
          -e $e -m $m \
          --cut-points-output ${path_cut_points_txt} | \
          pigz -c \
          > ${path_bed_gz}
        
        grep chm13 /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.fai | cut -f 1 | while read ref; do
          echo $ref
          \time -v ~/tools/odgi/bin/odgi-358537678174adc975415a488b458725d8a213be untangle -t 12 -P \
              -i $path_input_og \
              -r $ref \
              -e $e -m $m \
              --cut-points-input ${path_cut_points_txt} \
              -j 0 -n 100 | pigz -c > /lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$ref.e$e.m$m.j0.n100.bed.gz &
        done
        sleep 1h 30m # wait seems not to work properly
        wait
    fi;
  done
done
```

Grounding and plotting:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/untangle/grounded

path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.100kbps.pq_contigs.union.s100k.l300k.p98.n97/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.og
prefix=$(basename $path_input_og .og)

# Grounding (applying filters)
for e in 50000; do
  for m in 1000; do
      for j in 0 0.8; do
        j_str=$(echo $j | sed 's/\.//g')
        (seq 1 10; seq 15 5 100) | while read n; do 
            echo "-e $e -m $m -j $j -n $n"
            
            grep chm13 /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.fai | cut -f 1 | while read ref; do
              echo -e "\t"$ref
              
              path_grounded_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.tsv.gz
              
              if [[ ! -s ${path_grounded_tsv_gz} ]]; then
                  ( echo query query.begin query.end target target.begin target.end jaccard strand self.coverage ref ref.begin ref.end | tr ' ' '\t'
                    join \
                      <(zcat /lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.bed.gz | awk '{ print $1"_"$2, $0 }' | tr ' ' '\t' | sort -k 1,1) \
                      <(zcat /lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$ref.e$e.m$m.n100.j0.bed.gz | awk -v j=$j -v n=$n '{ if($7 >= j && $10 <= n) {print $1"_"$2, $0} }' | tr ' ' '\t' | sort -k 1,1) | \
                    tr ' ' '\t' | grep -v '^#' | cut -f 2- | cut -f -9,14-16 ) | tr ' ' '\t' | pigz -c > ${path_grounded_tsv_gz}
              fi;
            done
        done
      done
    done
done

# Take only reliable blocks (flagged with "Hh" or "Hc", https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/README.md#components)
#for e in 5000; do
#  for m in 1000; do
#      for j in 0.8; do
#        j_str=$(echo $j | sed 's/\.//g')
#        (seq 1 1) | while read n; do 
#            echo "-e $e -m $m -j $j -n $n"
#    
#            grep chm13 /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.fai | cut -f 1 | while read ref; do
#              echo -e "\t"$ref
#              
#              path_grounded_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.tsv.gz
#  
#              echo $path_grounded_tsv_gz
#             zcat $path_grounded_tsv_gz | awk '{print $0 >> $1".bed"}' example.bed
#            done
#        done
#      done
#    done
#done

# Take pq-untangling contigs
for e in 5000 10000 20000 30000 40000 50000 100000; do
  for m in 1000 10000; do
      for j in 0 0.8; do
        j_str=$(echo $j | sed 's/\.//g')
        (seq 1 10; seq 15 5 100) | while read n; do 
            echo "-e $e -m $m -j $j -n $n"
    
            grep chm13 /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.fai | cut -f 1 | while read ref; do
              echo $ref
              
              path_grounded_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.tsv.gz
              path_grounded_pq_touching_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv
              
              if [[ ! -s ${path_grounded_pq_touching_tsv}.gz ]]; then
                  # "p-touching contigs" intersected "q-touching contigs"
                  comm -12 \
                    <(bedtools intersect \
                      -a <(zcat ${path_grounded_tsv_gz} | awk -v OFS="\t" '{print $10,$11,$12,$1, "", "+"}' | grep query -v | bedtools sort) \
                      -b <(cat /lizardfs/guarracino/chromosome_communities/pq_contigs/p_arms.bed | bedtools sort) | \
                      #awk '$3-$2+1>=100000' | \
                      cut -f 4 | sort | uniq) \
                    <(bedtools intersect \
                      -a <(zcat ${path_grounded_tsv_gz} | awk -v OFS="\t" '{print $10,$11,$12,$1, "", "+"}' | grep query -v | bedtools sort) \
                      -b <(cat /lizardfs/guarracino/chromosome_communities/pq_contigs/q_arms.bed | bedtools sort) | \
                      #awk '$3-$2+1>=100000' | \
                      cut -f 4 | sort | uniq) | \
                      grep 'chm13#\|grch38#' -v > $ref.tmp.txt
                
                  # Put header (with a new 'target' column), take intersection, and re-add other acrocentric references
                  cat \
                    <(zcat ${path_grounded_tsv_gz} | head -n 1 | awk -v OFS='\t' '{print $0, "grounded.target"}') \
                    <(zgrep ${path_grounded_tsv_gz} -f $ref.tmp.txt | awk -v OFS='\t' -v ref=$ref '{print $0, ref}') \
                    <(zgrep '^grch38\|^chm13' ${path_grounded_tsv_gz}  | awk -v OFS='\t' -v ref=$ref '{print $0, ref}') \
                     > ${path_grounded_pq_touching_tsv}
                
                  # Add annotation
                  zgrep rDNA /lizardfs/guarracino/chromosome_communities/data/chm13.CentromeresAndTelomeres.CenSatAnnotation.txt.gz | \
                    sed 's/chr/chm13#chr/g' | \
                    grep $ref | \
                    awk -v OFS='\t' -v ref=$ref '{print $1"#"$4,".",".",$1,".",".",".",".",".",".",$2,$3, ref}' >> ${path_grounded_pq_touching_tsv}
                  grep acen /lizardfs/guarracino/chromosome_communities/data/chm13.CytoBandIdeo.v2.txt | \
                    bedtools merge | \
                    grep 'chr13\|chr14\|chr15\|chr21\|chr22' | \
                    sed 's/chr/chm13#chr/g'  | \
                    grep $ref | \
                    awk -v OFS='\t' -v ref=$ref '{print $1"#centromere",".",".",$1,".",".",".",".",".",".",$2,$3, ref}' >> ${path_grounded_pq_touching_tsv}
                  
                  pigz ${path_grounded_pq_touching_tsv}
              fi;
            done
        done
      done
    done
done
```

Plotting:

```shell
# Dependencies
guix install r
guix install r-ggplot2
guix install r-ggforce
guix install poppler # For pdfunite
```

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/
```

Single plots:

```shell
for e in 5000 10000 20000 30000 40000 50000 100000; do
  for m in 1000 10000; do
      for j in 0 0.8; do
        j_str=$(echo $j | sed 's/\.//g')
        (seq 1 10; seq 15 5 100) | while read n; do 
            echo "-e $e -m $m -j $j -n $n"
    
            grep chm13 /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.fai | cut -f 1 | while read ref; do
              echo $ref
              
              path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz
            
              Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle.R ${path_grounded_pq_touching_tsv_gz} "$ref -e $e -m $m -j $j -n $n"
            done
        done
      done
    done
done
```

Merged plots:

```shell
for e in 5000 10000 20000 30000 40000 50000 100000; do
  for m in 1000 10000; do
      for j in 0 0.8; do
        j_str=$(echo $j | sed 's/\.//g')
        (seq 1 10; seq 15 5 100) | while read n; do 
            echo "-e $e -m $m -j $j -n $n"
    
            path_grounded_pq_touching_all_chromosomes_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chrACRO.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz
            if [[ ! -s ${path_grounded_pq_touching_all_chromosomes_tsv_gz} ]]; then
                # Merge all acros together
                cat \
                  <(zcat /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.*.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz | head -n 1) \
                  <(zcat /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.*.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz | grep query -v) |\
                  pigz -c > x.tsv.gz
                mv x.tsv.gz ${path_grounded_pq_touching_all_chromosomes_tsv_gz}
            fi;

            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_all.R ${path_grounded_pq_touching_all_chromosomes_tsv_gz} "-e $e -m $m -j $j -n $n"
        done
      done
    done
done

mv /lizardfs/guarracino/chromosome_communities/untangle/grounded/*.pdf /lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/

for e in 5000 10000 20000 30000 40000 50000 100000; do
  for m in 1000 10000; do
      for j in 0 0.8; do
        echo "-e $e -m $m -j $j"
        
        j_str=$(echo $j | sed 's/\.//g')
        PDFs=$((seq 1 10; seq 15 5 100) | while read n; do \
          echo /lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/$prefix.untangle.chm13#chrACRO.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz.pdf
        done | tr '\n' ' ')
        #echo $PDFs
        
        # Too slow
        #gs -sDEVICE=pdfwrite \
        #  -dNOPAUSE -dBATCH -dSAFER \
        #  -sOutputFile=/lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/$prefix.untangle.chm13#chrACRO.e$e.m$m.j${j_str}.grounded.pq_touching.tsv.gz.pdf \
        #  $PDFs
          
        /gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite $PDFs /lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/$prefix.untangle.chm13#chrACRO.e$e.m$m.j${j_str}.merged.grounded.pq_touching.tsv.gz.pdf
      done
    done
done
```





```shell
( echo A | tr ' ' '\n') | while read i; do
  sbatch -p headnode -w octopus01 -c 32 --wrap 'hostname; cd /scratch && /gnu/store/swnkjnc9wj6i1cl9iqa79chnf40r1327-pggb-0.2.0+640bf6b-5/bin/pggb -i /lizardfs/erikg/HPRC/year1v2genbank/parts/chr'$i'.pan.fa -o chr'$i'.pan.s100k.l300k.p97.n200 -t 32 -p 97 -s 100000 -l 300000 -n 200 -k 311 -G 13117,13219 -O 0.03 -T 32 -v -V chm13:#,grch38:# --resume; mv /scratch/chr'$i'.pan /lizardfs/guarracino/chromosome_communities/graphs';
done # >>pggb.jobids
```






Prepare the references:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/recombination
cd /lizardfs/guarracino/chromosome_communities/recombination

samtools faidx /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa \
  $(grep 'chr13\|chr14\|chr15\|chr21\|chr22' /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa.fai | grep '_' -v | cut -f 1) > chm13.ACRO.fa
samtools faidx /lizardfs/erikg/HPRC/year1v2genbank/assemblies/grch38.fa \
  $(grep 'chr13\|chr14\|chr15\|chr21\|chr22' /lizardfs/erikg/HPRC/year1v2genbank/assemblies/grch38.fa.fai | grep '_' -v | cut -f 1) > grch38.ACRO.fa
```







[comment]: <> (```shell)

[comment]: <> (perl ~/git/vcf-conversion-tools/vcf2MS.pl example.vcf example.ms 629)

[comment]: <> (perl ~/git/vcf-conversion-tools/MS2LDhat.pl example.ms example.ldhat 629)

[comment]: <> (#or)

[comment]: <> (perl ~/git/vcf-conversion-tools/vcf2MS.pl example.10samples.vcf example.ms 11)

[comment]: <> (perl ~/git/vcf-conversion-tools/MS2LDhat.pl example.ms example 11)

[comment]: <> (~/git/LDhat/pairwise -seq example.ldhat.sites -loc example.ldhat.locs)

[comment]: <> (pip3 install biopython)

[comment]: <> (pip3 install pyfaidx)

[comment]: <> (pip3 install PyVCF)

[comment]: <> (python3 vcf2alignedFasta.py example.10samples.vcf.gz chr17.fasta )

[comment]: <> (# hack output &#40;num seq, length seq, 1 &#40;haploid&#41; / 2 &#40;diploid&#41;)

[comment]: <> (~/git/LDhat/convert -seq example.fasta )

[comment]: <> (~/git/LDhat/pairwise -seq sites.txt -loc locs.txt )

[comment]: <> (path_vcf=example.10samples.vcf)

[comment]: <> (path_vcf=pggb.wgg.88.chm13.1-22+X.norm.max50.vcf)

[comment]: <> (vk phylo fasta ${path_vcf} > ${path_vcf}.tmp.fasta # todo to improve)

[comment]: <> (seq_num=$&#40;grep '^>' ${path_vcf}.tmp.fasta -c&#41;)

[comment]: <> (seq_length=$&#40;head ${path_vcf}.tmp.fasta -n 2 | tail -n 1 | awk '{print length&#40;$0&#41;}'&#41;)

[comment]: <> (# 1 &#40;haploid&#41; / 2 &#40;diploid&#41;)

[comment]: <> (cat <&#40;echo -e $seq_num"\t"$seq_length"\t"1&#41; <&#40;cat ${path_vcf}.tmp.fasta&#41; > ${path_vcf}.fasta)

[comment]: <> (rm ${path_vcf}.tmp.fasta)

[comment]: <> (~/git/LDhat/convert -seq ${path_vcf}.fasta)

[comment]: <> (~/git/LDhat/pairwise -seq sites.txt -loc locs.txt )

[comment]: <> (```)

# OLD STUFF (TO DELETE OR FINISH)

Prepare CHM13's acrocentric chromosomes:

```shell
samtools faidx chm13.fa.gz $(grep 'chr13\|chr14\|chr15\|chr21\|chr22' chm13.fa.gz.fai | cut -f 1) | bgzip -@ 48 -c > chm13.ACRO.fa.gz && samtools faidx chm13.ACRO.fa.gz
rm chm13.draft_v1.1.fasta.gz chm13.fa.gz*
```

## Pangenome building

# todo: to take the last chrY (check the differences)

Use HG002's chrY as an alternative reference, as GRCh38's chrY is incomplete. Include also HG002's chrX.

```shell
mkdir -p /lizardfs/guarracino/HPRC/chromosome_communities/data
cd /lizardfs/guarracino/HPRC/chromosome_communities/data

#2021-08-27T19:04:55.000Z        3.0 GB         v2.fasta
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/T2T/HG002XY/v2/v2.fasta
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/T2T/HG002XY/v2/v2.fasta.fai
samtools faidx v2.fasta $(grep hg002 v2.fasta.fai | cut -f 1) | \
 sed 's/>chrX_hg002_v2/>HG002#chrX/g' | sed 's/>chrY_hg002_v2/>HG002#chrY/g' > H002.chrXY.fa

rm v2.fasta*
```

Put HG002's chrX and chrY with the partitioned chrXs and chrYs.

```shell
# Prepare sequence order, with all references on the top
grep chr /lizardfs/erikg/HPRC/year1v2genbank/parts/chrS.pan.fa.fai | cut -f 1 > sequence_order.txt
grep chr /lizardfs/erikg/HPRC/year1v2genbank/parts/chrS.pan.fa.fai -v | cut -f 1 >> sequence_order.txt
rm sequence_order.txt

cat H002.chrXY.fa <(samtools faidx /lizardfs/erikg/HPRC/year1v2genbank/parts/chrS.pan.fa $(cat sequence_order.txt)) | bgzip -@ 48 > chrS.pan+HG002chrXY.fa.gz
samtools faidx chrS.pan+HG002chrXY.fa.gz

cd /lizardfs/guarracino/HPRC/chromosome_communities/
```

### Build graphs

Apply `pggb` on the chromosome-partitioned HPRC dataset. We make two graphs considering the following chromosome
jointly:

- sex chromosomes (chrX and chrY).

```shell
sbatch -p workers -c 48 -w octopus11 --wrap 'hostname; cd /scratch && /gnu/store/swnkjnc9wj6i1cl9iqa79chnf40r1327-pggb-0.2.0+640bf6b-5/bin/pggb -i /lizardfs/guarracino/HPRC/chromosome_communities/data/chrS.pan+HG002chrXY.fa.gz -o chrS.pan+HG002chrXY -t 48 -p 98 -s 100000 -n 90 -k 311 -O 0.03 -T 48 -U -v -L -V chm13:#,grch38:# -Z ; mv /scratch/chrS.pan+HG002chrXY /lizardfs/guarracino/HPRC/chromosome_communities/'
```

- acrocentric chromosomes (chr13, chr14, chr15, chr21, and chr22);

```shell
( echo A | tr ' ' '\n') | while read i; do sbatch -p workers -c 32 -w octopus02 --wrap 'hostname; cd /scratch && /gnu/store/swnkjnc9wj6i1cl9iqa79chnf40r1327-pggb-0.2.0+640bf6b-5/bin/pggb -i /lizardfs/erikg/HPRC/year1v2genbank/parts/chr'$i'.pan.fa -o chr'$i'.pan.s100k.l300k.p97.n200 -t 32 -p 97 -s 100000 -l 300000 -n 200 -k 311 -O 0.03 -T 32 -U -v -L -V chm13:#,grch38:#; mv /scratch/chr'$i'.pan /lizardfs/guarracino/chromosome_communities/graphs'; done # >>pggb.jobids
```

### Untangling

Untangle all contigs in the sex graph by using the GRCh38 reference and the new HG002 _de novo_ assembly:

```shell
bash untangle.sh chrS.pan+HG002chrY/chrS.pan+HG002chrY.fa.gz.a2fb268.4030258.6a1ecc2.smooth.og.gz grch38 "chrX chrY" 10000 sex
bash untangle.sh chrS.pan+HG002chrY/chrS.pan+HG002chrY.fa.gz.a2fb268.4030258.6a1ecc2.smooth.og.gz "chm13 HG002" "chrX chrY" 10000 sex
```

### ...

ToDo