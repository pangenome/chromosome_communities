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


## Collect contigs running from the p-arm to the q-arm of the acrocentric chromosomes

Prepare CHM13's acrocentric chromosomes:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/assemblies
cd /lizardfs/guarracino/chromosome_communities/assemblies

wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/CHM13/assemblies/chm13.draft_v1.1.fasta.gz
( ~/tools/fastix/target/release/fastix-331c1159ea16625ee79d1a82522e800c99206834 -p 'chm13#' <(zcat chm13.draft_v1.1.fasta.gz) | bgzip -@ 48 -c >chm13.fa.gz && samtools faidx chm13.fa.gz)
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
        sbatch -p allnodes -c 24 --job-name acrovsref --wrap 'hostname; cd /scratch; \time -v /gnu/store/d6ajy0ajxkbb22gkgi32g5c4jy8rs8bv-wfmash-0.6.0+948f168-21/bin/wfmash /lizardfs/guarracino/chromosome_communities/assemblies/chm13.chr'$f'.fa.gz /lizardfs/erikg/HPRC/year1v2genbank/parts/chr'$f'.pan.fa -s '$s' -l '$l' -p '$p' -n 1 -t 24 -m > /lizardfs/guarracino/chromosome_communities/mappings/chr'$f'.vs.chm13.chr'$f'.s'$s'.l'$l'.p'$p'.n1.paf'
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
   sort | bedtools complement -i - -g <(cut -f 1,2 /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa.fai | sort) | \
   grep 'chr13\|chr14\|chr15\|chr21\|chr22' | \
   sed -n 1~2p > /lizardfs/guarracino/chromosome_communities/pq_contigs/p_arms.bed

bedtools slop \
  -i <(grep acen /lizardfs/guarracino/chromosome_communities/data/chm13.CytoBandIdeo.v2.txt | bedtools merge | sed 's/chr/chm13#chr/g' | sort) \
  -g <(cut -f 1,2 /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa.fai | sort) \
  -b 1000000 | \
   sort | bedtools complement -i - -g <(cut -f 1,2 /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa.fai | sort) | \
   grep 'chr13\|chr14\|chr15\|chr21\|chr22' | \
   sed -n 2~2p > /lizardfs/guarracino/chromosome_communities/pq_contigs/q_arms.bed
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

# Debugging
#for s in 300k 200k 150k 100k 50k 20k ; do
#  for p in 98 95 90 85 80; do
#    s_no_k=${s::-1}
#    l_no_k=$(echo $s_no_k '*' 3 | bc)
#    l=${l_no_k}k
#    
#    (seq 13 15; seq 21 22) | while read f; do      
#      if [[ ! -s /lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.chr$f.s$s.l$l.p$p.n1.pq_contigs.txt ]]; then
#          echo s$s l$l p$p chr$f
#          
#          bedtools intersect \
#              -a <(cat /lizardfs/guarracino/chromosome_communities/mappings/chr$f.vs.chm13.chr$f.s$s.l$l.p$p.n1.paf | awk -v OFS='\t' '{print $6, $8, $9, $1, "", "+"}' | bedtools sort) \
#              -b <(grep chr$f /lizardfs/guarracino/chromosome_communities/pq_contigs/p_arms.bed | bedtools sort) | \
#              awk '$3-$2+1>=100000' > /lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.chr$f.s$s.l$l.p$p.n1.p_contigs.100kbps.bed
#          bedtools intersect \
#              -a <(cat /lizardfs/guarracino/chromosome_communities/mappings/chr$f.vs.chm13.chr$f.s$s.l$l.p$p.n1.paf | awk -v OFS='\t' '{print $6, $8, $9, $1, "", "+"}' | bedtools sort) \
#              -b <(grep chr$f /lizardfs/guarracino/chromosome_communities/pq_contigs/q_arms.bed | bedtools sort) | \
#              awk '$3-$2+1>=100000' > /lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.chr$f.s$s.l$l.p$p.n1.q_contigs.100kbps.bed
#              
#          comm -12 \
#            <(cut -f 4 /lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.chr$f.s$s.l$l.p$p.n1.p_contigs.100kbps.bed | sort | uniq) \
#            <(cut -f 4 /lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.chr$f.s$s.l$l.p$p.n1.q_contigs.100kbps.bed | sort | uniq) | \
#              grep 'chm13#\|grch38#' -v > /lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.chr$f.s$s.l$l.p$p.n1.pq_contigs.100kbps.txt
#              
#              
#          cat /lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.chr$f.s$s.l$l.p$p.n1.pq_contigs.100kbps.txt | while read contig; do
#            leftmost_coordinate=$(grep $contig /lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.chr$f.s$s.l$l.p$p.n1.q_contigs.100kbps.bed | cut -f 2,3 | tr '\t' '\n' | sort -k 1n | head -n 1)
#            rightmost_coordinate=$(grep $contig /lizardfs/guarracino/chromosome_communities/pq_contigs/chr$f.vs.chm13.chr$f.s$s.l$l.p$p.n1.q_contigs.100kbps.bed | cut -f 2,3 | tr '\t' '\n' | sort -k 1n | tail -n 1)
#            
#            mapping_length_on_the_reference=$(echo $rightmost_coordinate-$leftmost_coordinate+1 | bc)
#            contig_length=$(grep $contig /lizardfs/erikg/HPRC/year1v2genbank/parts/chr$f.pan.fa.fai | cut -f 2)
#            
#            ratio=$(echo "scale=4; $mapping_length_on_the_reference / $contig_length" | bc)
#            echo $contig $mapping_length_on_the_reference $contig_length $ratio | tr ' ' '\t'
#          done | column -t
#      fi;
#    done
#  done
#done
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

Put all together with the reference:

```shell
path_grch38_acro_fa_gz=/lizardfs/guarracino/chromosome_communities/pq_contigs/grch38.ACRO.fa.gz
samtools faidx /lizardfs/erikg/HPRC/year1v2genbank/assemblies/grch38.fa $(grep 'chr13\|chr14\|chr15\|chr21\|chr22' /lizardfs/erikg/HPRC/year1v2genbank/assemblies/grch38.fa.fai | grep '_' -v | cut -f 1) | bgzip -@ 48 -c > $path_grch38_acro_fa_gz

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
mkdir -p /lizardfs/guarracino/chromosome_communities/untangle

path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.100kbps.pq_contigs.union.s100k.l300k.p98.n97/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.og
prefix=$(basename $path_input_og .og)

/gnu/store/2ln6zv8mk6aqqzcib39rgi11x2hn7mv9-odgi-0.6.2+9e9c481-13/bin/odgi untangle -t 32 -P \
  -i $path_input_og \
  -R <(grep chm13 /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.fai | cut -f 1) \
  -e 50000 -m 1000 > /lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e50000.m1000.tsv

grep chm13 /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.fai | cut -f 1 | while read ref; do
  echo $ref
  /gnu/store/2ln6zv8mk6aqqzcib39rgi11x2hn7mv9-odgi-0.6.2+9e9c481-13/bin/odgi untangle -t 32 -P \
      -i $path_input_og \
      -r $ref \
      -e 50000 -m 1000 \
      -n 50 -j 0.8 > /lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$ref.e50000.m1000.n50.j08.tsv
done

grep chm13 /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.fai | cut -f 1 | while read ref; do
  echo $ref
    ( echo query query.begin query.end target target.begin target.end jaccard strand self.coverage ref ref.begin ref.end | tr ' ' '\t'
    join \
      <(cat /lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e50000.m1000.tsv | awk '{ print $1"_"$2, $0 }' | tr ' ' '\t' | sort -k 1,1) \
      <(cat /lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$ref.e50000.m1000.n50.j08.tsv | awk '{ print $1"_"$2, $0 }' | tr ' ' '\t' | sort -k 1,1) | \
      tr ' ' '\t' | grep -v '^#' | cut -f 2- | cut -f -9,14-16 ) | tr ' ' '\t' > /lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$ref.e50000.m1000.n50.j08.grounded.tsv
done
```

Plotting:

```R
library(ggplot2)
x <- read.delim('/home/guarracino/Downloads/Pangenomics/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.untangle.chm13#chr22.e5000.m1000.n50.j08.grounded.tsv')
p <- ggplot(x, aes(x=ref.begin+(ref.end-ref.begin)/2, width=ref.end-ref.begin-200, y=query, fill=target)) + geom_tile()
ggsave(plot = p, 'chr22.e5000.m1000.n50.j08.grounded.pdf', width = 21, height = 14,  units = "cm", dpi = 300,  bg = "transparent")
```


Take the contigs that go from one side of the centromere to the other:

```shell

```

### Data preparation

ls HG002's chrY as an alternative reference, as GRCh38's chrY is incomplete. Include also HG002's chrX.

```
mkdir -p /lizardfs/guarracino/HPRC/chromosome_communities/data
cd /lizardfs/guarracino/HPRC/chromosome_communities/data

wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/T2T/HG002XY/v2/v2.fasta
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/T2T/HG002XY/v2/v2.fasta.fai
samtools faidx v2.fasta $(grep hg002 v2.fasta.fai | cut -f 1) | \
 sed 's/>chrX_hg002_v2/>HG002#chrX/g' | sed 's/>chrY_hg002_v2/>HG002#chrY/g' > H002.chrXY.fa

rm v2.fasta*
```

Put HG002's chrX and chrY with the partitioned chrXs and chrYs.

```
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

```
sbatch -p workers -c 48 -w octopus11 --wrap 'hostname; cd /scratch && /gnu/store/2rwrch6gc5r5pazikhfc25j2am2rh22a-pggb-0.2.0+531f85f-1/bin/pggb -i /lizardfs/guarracino/HPRC/chromosome_communities/data/chrS.pan+HG002chrXY.fa.gz -o chrS.pan+HG002chrXY -t 48 -p 98 -s 100000 -n 90 -k 311 -O 0.03 -T 48 -U -v -L -V chm13:#,grch38:# -Z ; mv /scratch/chrS.pan+HG002chrXY /lizardfs/guarracino/HPRC/chromosome_communities/'
```

- acrocentric chromosomes (chr13, chr14, chr15, chr21, and chr22);

```
( echo A | tr ' ' '\n') | while read i; do sbatch -p workers -c 48 -w octopus02 --wrap 'hostname; cd /scratch && /gnu/store/2rwrch6gc5r5pazikhfc25j2am2rh22a-pggb-0.2.0+531f85f-1/bin/pggb -i /lizardfs/erikg/HPRC/year1v2genbank/parts/chr'$i'.pan.fa -o chr'$i'.pan -t 48 -p 98 -s 100000 -n 90 -k 311 -O 0.03 -T 48 -U -v -L -V chm13:#,grch38:# -Z ; mv /scratch/chr'$i'.pan /lizardfs/guarracino/HPRC/chromosome_communities/'; done # >>pggb.jobids
```

### Untangling

Untangle all contigs in the sex graph by using the GRCh38 reference and the new HG002 _de novo_ assembly:

```
bash untangle.sh chrS.pan+HG002chrY/chrS.pan+HG002chrY.fa.gz.a2fb268.4030258.6a1ecc2.smooth.og.gz grch38 "chrX chrY" 10000 sex
bash untangle.sh chrS.pan+HG002chrY/chrS.pan+HG002chrY.fa.gz.a2fb268.4030258.6a1ecc2.smooth.og.gz "chm13 HG002" "chrX chrY" 10000 sex
```


### ...

ToDo
