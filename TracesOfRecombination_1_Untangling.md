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
git checkout 694948ccf31e7b565449cc056668e9dcc8cc0a3e
cmake -H. -Bbuild && cmake --build build -- -j 48
mv bin/odgi bin/odgi-694948ccf31e7b565449cc056668e9dcc8cc0a3e
cd ..

(echo pggb wfmash seqwish smoothxg odgi gfaffix | tr ' ' '\n') | while read tool; do ls -l $(which $tool); done | cut -f 13 -d ' '
/gnu/store/swnkjnc9wj6i1cl9iqa79chnf40r1327-pggb-0.2.0+640bf6b-5/bin/pggb
/gnu/store/d6ajy0ajxkbb22gkgi32g5c4jy8rs8bv-wfmash-0.6.0+948f168-21/bin/wfmash
/gnu/store/v1q6ja2fy3c7fcz7bnr6k3x54amynhyp-seqwish-0.7.3+51ee550-1/bin/seqwish
/gnu/store/p8wflvwwxiih4z9ds5hfkin8mjld6qw2-smoothxg-0.6.0+0f383e5-10/bin/smoothxg
/gnu/store/2ln6zv8mk6aqqzcib39rgi11x2hn7mv9-odgi-0.6.2+9e9c481-13/bin/odgi
/gnu/store/ccw48k7h8v1brz3ap0sj3bcwvvmk6xra-gfaffix-0.1.2.2/bin/gfaffix

vg
vg: variation graph tool, version v1.38.0 "Canossa"
```

## Preparation

Clone the repository:

```shell
cd /lizardfs/guarracino/
git clone --recursive https://github.com/pangenome/chromosome_communities.git
```

Prepare verkko's HG002 contigs (HiFi+ONT-based):

```shell
cd /lizardfs/guarracino/chromosome_communities/assemblies

### Get the "hg002-prox.fna" file from Globus (t2t-share/hg002/acros/beta folder)

# Rename contigs
# Ugly! ($f doesn't work)
cat hg002-prox.fna | \
  sed "s/chr13_prox_MAT/HG002#MAT#chr13.prox/g" |\
  sed "s/chr13_prox_PAT/HG002#PAT#chr13.prox/g" |\
  sed "s/chr14_prox_MAT/HG002#MAT#chr14.prox/g" |\
  sed "s/chr14_prox_PAT/HG002#PAT#chr14.prox/g" |\
  sed "s/chr15_prox_MAT/HG002#MAT#chr15.prox/g" |\
  sed "s/chr15_prox_PAT/HG002#PAT#chr15.prox/g" |\
  sed "s/chr21_prox_MAT/HG002#MAT#chr21.prox/g" |\
  sed "s/chr21_prox_PAT/HG002#PAT#chr21.prox/g" |\
  sed "s/chr22_prox_MAT/HG002#MAT#chr22.prox/g" |\
  sed "s/chr22_prox_PAT/HG002#PAT#chr22.prox/g" > hg002-prox.renamed.fna
samtools faidx hg002-prox.renamed.fna

# Partition by chromosomes
(seq 13 15; seq 21 22) | while read f; do
  samtools faidx hg002-prox.renamed.fna $(grep chr$f hg002-prox.renamed.fna.fai | cut -f 1) |\
    bgzip -@ 48 -c > hg002.chr$f.prox.fa.gz
  samtools faidx hg002.chr$f.prox.fa.gz
done
```

Prepare HiFi-based HG002 contigs:

```shell
# Download
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/assemblies/year1_freeze_assembly_v2.1/HG002.maternal.f1_assembly_v2.1.fa.gz
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/assemblies/year1_freeze_assembly_v2.1/HG002.paternal.f1_assembly_v2.1.fa.gz

for hap in maternal paternal; do
  echo $hap
  gunzip HG002.$hap.f1_assembly_v2.1.fa.gz -c | bgzip -@ 48 -c > HG002.$hap.fa.gz
  samtools faidx HG002.$hap.fa.gz
done
rm HG002.*.f1_assembly_v2.1.fa.gz

# HG002 has low quality, so no pq-contigs. We use verkko's assemblies to collect all its acro-centric contigs
mkdir -p partitioning

for hap in mat pat; do
  samtools faidx hg002-prox.renamed.fna $(grep $hap hg002-prox.renamed.fna.fai -i | cut -f 1) | bgzip -@ 48 -c > hg002-prox.renamed.$hap.fna.gz
  samtools faidx hg002-prox.renamed.$hap.fna.gz
done

## Map against the verkko's assemblies (more strict identity threshold)
for hap in mat pat; do
  reffa=/lizardfs/guarracino/chromosome_communities/assemblies/hg002-prox.renamed.$hap.fna.gz
  hapfa=/lizardfs/guarracino/chromosome_communities/assemblies/HG002.${hap}ernal.fa.gz 
  
  sbatch -c 24 -p workers --job-name HG002 --wrap 'hostname; cd /scratch; \time -v /gnu/store/d6ajy0ajxkbb22gkgi32g5c4jy8rs8bv-wfmash-0.6.0+948f168-21/bin/wfmash -t 24 -m -N -s 50k -p 95 '$reffa' '$hapfa' > /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/HG002.'$hap'.vs.hg002-prox.'$hap'.paf'
done

## Collect unmapped contigs and remap them in split mode
for hap in mat pat; do
  echo $hap
  hapfa=/lizardfs/guarracino/chromosome_communities/assemblies/HG002.${hap}ernal.fa.gz
  nonsplitpaf=/lizardfs/guarracino/chromosome_communities/assemblies/partitioning/HG002.$hap.vs.hg002-prox.$hap.paf

  comm -23 <(cut -f 1 $hapfa.fai | sort) <(cut -f 1 $nonsplitpaf | sort) > partitioning/HG002.$hap.vs.hg002-prox.$hap.missing.txt
  if [[ $(wc -l partitioning/HG002.$hap.vs.hg002-prox.$hap.missing.txt | cut -f 1 -d\ ) != 0 ]];
  then
    missingfa=/lizardfs/guarracino/chromosome_communities/assemblies/partitioning/HG002.$hap.vs.hg002-prox.$hap.missing.fa.gz
    
    samtools faidx $hapfa $(tr '\n' ' ' < partitioning/HG002.$hap.vs.hg002-prox.$hap.missing.txt) | bgzip -@ 48 -c > $missingfa
    samtools faidx $missingfa
    
    sbatch -c 24 -p workers --job-name HG002 --wrap 'hostname; cd /scratch; \time -v /gnu/store/d6ajy0ajxkbb22gkgi32g5c4jy8rs8bv-wfmash-0.6.0+948f168-21/bin/wfmash -t 24 -m -s 50k -p 95 '$reffa' '$missingfa' > /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/HG002.'$hap'.split.vs.hg002-prox.'$hap'.paf'
  fi
done

# Nothing recoverable
cat partitioning/HG002.*.split.vs.hg002-prox.*.paf

## Subset by acrocentric chromosome
( seq 13 15; seq 21 22 ) | while read f; do
  echo $f
  
  for hap in mat pat; do
    hapfa=/lizardfs/guarracino/chromosome_communities/assemblies/HG002.${hap}ernal.fa.gz
  
    samtools faidx $hapfa $(awk '$6 ~ "chr'$f'.prox$"' partitioning/HG002.$hap.vs.hg002-prox.$hap.paf | cut -f 1) >> hg002.chr$f.hifi.fa
  done
  
  bgzip -@ 48 hg002.chr$f.hifi.fa
  samtools faidx hg002.chr$f.hifi.fa.gz
done
```

[comment]: <> (## Get coverage analysis BED files)
[comment]: <> (```shell)
[comment]: <> (mkdir -p /lizardfs/guarracino/chromosome_communities/coverage_analysis_y1_genbank)
[comment]: <> (cd /lizardfs/guarracino/chromosome_communities/coverage_analysis_y1_genbank)
[comment]: <> (prefix=https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/e9ad8022-1b30-11ec-ab04-0a13c5208311--COVERAGE_ANALYSIS_Y1_GENBANK/FLAGGER/JAN_09_2022/FINAL_HIFI_BASED/FLAGGER_HIFI_ASM_BEDS/)
[comment]: <> (cut -f 1 /lizardfs/erikg/HPRC/year1v2genbank/parts/HPRCy1.pan.fa.gz.fai | cut -f 1 -d '#' | grep 'chm13\|grch38' -v | sort | uniq | while read sample; do)
[comment]: <> (  wget -c $prefix$sample.hifi.flagger_final.bed)
[comment]: <> (done)
[comment]: <> (```)

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

Prepare the BED files for p-arms and q-arms (-/+ 1Mbps):

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/pq_contigs
cd /lizardfs/guarracino/chromosome_communities/pq_contigs

bedtools slop \
    -i <(sed 's/chr/chm13#chr/g' /lizardfs/guarracino/chromosome_communities/data/chm13.centromeres.approximate.bed | bedtools sort) \
    -g <(cut -f 1,2 /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa.fai | sort) \
    -b 1000000 | \
  bedtools sort | \
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

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/graphs

num_of_haplotypes=$(cat /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.fa.gz.fai | wc -l)
sbatch -p highmem -c 48 --job-name acropggb --wrap 'hostname; cd /scratch && /gnu/store/swnkjnc9wj6i1cl9iqa79chnf40r1327-pggb-0.2.0+640bf6b-5/bin/pggb -i /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.fa.gz -o chrACRO+refs.1kbps.pq_contigs.union.hg002prox.s100k.l300k.p98.n'$num_of_haplotypes' -t 48 -s 100k -l 300k -p 98 -n '$num_of_haplotypes' -k 311 -G 13117,13219 -O 0.03 -T 48 -v -V chm13:#,grch38:#; mv /scratch/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.s100k.l300k.p98.n'$num_of_haplotypes' /lizardfs/guarracino/chromosome_communities/graphs';

num_of_haplotypes=$(cat /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.fa.gz.fai | wc -l)
num_of_haplotypes_plus_a_bit=$(echo "$num_of_haplotypes + 10" | bc) # 5 mat acros + 5 pat acros
sbatch -p highmem -c 48 --job-name acropggb --wrap 'hostname; cd /scratch && /gnu/store/swnkjnc9wj6i1cl9iqa79chnf40r1327-pggb-0.2.0+640bf6b-5/bin/pggb -i /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.hg002hifi.fa.gz -o chrACRO+refs.1kbps.pq_contigs.union.hg002prox.hg002hifi.s100k.l300k.p98.n'$num_of_haplotypes_plus_a_bit' -t 48 -s 100k -l 300k -p 98 -n '$num_of_haplotypes_plus_a_bit' -k 311 -G 13117,13219 -O 0.03 -T 48 -v -V chm13:#,grch38:#; mv /scratch/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.hg002hifi.s100k.l300k.p98.n'$num_of_haplotypes_plus_a_bit' /lizardfs/guarracino/chromosome_communities/graphs';


#num_of_haplotypes=$(cut -f 1,2 -d '#' /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.fai | sort | uniq | wc -l)
#sbatch -p highmem -c 48 --job-name acropggb --wrap 'hostname; cd /scratch && /gnu/store/swnkjnc9wj6i1cl9iqa79chnf40r1327-pggb-0.2.0+640bf6b-5/bin/pggb -i /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz -o chrACRO+refs.50kbps.pq_contigs.union.hg002prox.s100k.l300k.p98.n'$num_of_haplotypes' -t 48 -s 100k -l 300k -p 98 -n '$num_of_haplotypes' -k 311 -G 13117,13219 -O 0.03 -T 48 -v -V chm13:#,grch38:#; mv /scratch/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.s100k.l300k.p98.n'$num_of_haplotypes' /lizardfs/guarracino/chromosome_communities/graphs';
```

### Untangling

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/untangle

path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.s100k.l300k.p98.n178/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.fa.gz.de31bcf.4030258.2385969.smooth.fix.og
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.hg002hifi.s100k.l300k.p98.n188/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.hg002hifi.fa.gz.4817bf7.4030258.57ed14c.smooth.fix.og

prefix=$(basename $path_input_og .og)
run_odgi=/home/guarracino/tools/odgi/bin/odgi-694948ccf31e7b565449cc056668e9dcc8cc0a3e

path_targets_txt=/lizardfs/guarracino/chromosome_communities/untangle/chm13.target_paths.txt
grep chm13 /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.fa.gz.fai | cut -f 1 > $path_targets_txt

# All references and emit cut points
for e in 5000 50000 100000; do
  for m in 500 1000 10000; do
    echo "-e $e -m $m"
    
    path_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.bed.gz
    path_cut_points_txt=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.cut_points.txt
    
    if [[ ! -s ${path_cut_points_txt} ]]; then
      sbatch -p workers -c 24 --job-name acrountangle --wrap '\time -v '$run_odgi' untangle -t 24 -P -i '$path_input_og' -R '$path_targets_txt' -e '$e' -m '$m' --cut-points-output '$path_cut_points_txt' -j 0 -n 1 | pigz -c > '$path_bed_gz';'
    fi;
  done
done

# Single reference by using the same cut points
for e in 5000 50000 100000; do
  for m in 500 1000 10000; do
    echo "-e $e -m $m"
    
    path_cut_points_txt=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.cut_points.txt
      
    cat $path_targets_txt | while read ref; do
      echo $ref
      
      path_ref_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$ref.e$e.m$m.j0.n100.bed.gz
      if [[ ! -s ${path_ref_bed_gz} ]]; then
        path_cut_points_txt=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.cut_points.txt
        
        sbatch -p workers -c 12 --job-name acrountangle --wrap '\time -v '$run_odgi' untangle -t 12 -P -i '$path_input_og' -r '$ref' -e '$e' -m '$m' --cut-points-input '$path_cut_points_txt' -j 0 -n 100 | pigz -c > '$path_ref_bed_gz';'
      fi;
    done
  done
done
```

Grounding (applying filters) and plotting:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/untangle/grounded

for e in 5000 50000 100000; do
  for m in 500 1000 10000; do
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

### Take only reliable blocks (flagged with "Hh" or "Hc", https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/README.md#components)???

# Take pq-untangling contigs
for e in 5000 50000 100000; do
  for m in 500 1000 10000; do
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

            comm -12 \
              <(sort $ref.tmp.txt) \
              <(cut -f 1 /lizardfs/guarracino/chromosome_communities/pq_contigs/`echo $ref | sed 's/chm13#//g'`.vs.chm13.1kbps.pq_contigs.union.hg002prox.fa.gz.fai | sort)  \
              > $ref.tmp2.txt
            rm $ref.tmp.txt
            
            # To plot the short (bad) HiFi-based HG002 contigs
            #cut -f 1 /lizardfs/guarracino/chromosome_communities/assemblies/hg002.`echo $ref | sed 's/chm13#//g'`.hifi.fa.gz.fai >> $ref.tmp2.txt         
               
            # Put header (with a new 'target' column), take intersection, and re-add other acrocentric references
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
              awk -v OFS='\t' -v ref=$ref '{print $1"#"$4,".",".",$1,".",".",".",".",".",".",$2,$3, ref}' >> ${path_grounded_pq_touching_tsv}
            cat /lizardfs/guarracino/chromosome_communities/data/chm13.centromeres.approximate.bed | \
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

Merged plots:

```shell
# Dependencies
#guix install r
#guix install r-ggplot2
#guix install r-ggforce
#guix install poppler # For pdfunite

mkdir -p /lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/

for e in 5000 50000 100000; do
  for m in 500 1000 10000; do
    for j in 0 0.8 0.95; do
        j_str=$(echo $j | sed 's/\.//g')
        (seq 1 5; seq 10 10 50) | while read n; do 
            echo "-e $e -m $m -j $j -n $n"
    
            path_grounded_pq_touching_all_chromosomes_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chrACRO.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz
            if [[ ! -s ${path_grounded_pq_touching_all_chromosomes_tsv_gz} ]]; then
                # Merge single reference results
                cat \
                  <(zcat /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.*.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz | head -n 1) \
                  <(zcat /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.*.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz | grep query -v) |\
                  pigz -c > x.tsv.gz
                # Rename after to avoid getting itself with the previous '*' expansion
                mv x.tsv.gz ${path_grounded_pq_touching_all_chromosomes_tsv_gz}
            fi;

            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_all.R ${path_grounded_pq_touching_all_chromosomes_tsv_gz} "-e $e -m $m -j $j -n $n" 25000000 200
        done
      done
    done
done

mv /lizardfs/guarracino/chromosome_communities/untangle/grounded/*.pdf /lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/

# Merge
for e in 5000 50000 100000; do
  for m in 500 1000 10000; do
    for j in 0 0.8 0.95; do
        echo "-e $e -m $m -j $j"
        
        j_str=$(echo $j | sed 's/\.//g')
        PDFs=$((seq 1 5; seq 10 10 50) | while read n; do \
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


## Sex chromosome

Use HG002's chrY as an alternative reference, as GRCh38's chrY is incomplete. Include also HG002's chrX.

```shell
cd /lizardfs/guarracino/chromosome_communities/assemblies

# Get HG002 chromosome X and Y
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/21edcb42-02c4-4e9f-b226-6773e62484a4--RU-HG002-commons/assembly/curated_round2/HG002.mat.cur.20211005.fasta.gz
gunzip HG002.mat.cur.20211005.fasta.gz
samtools faidx HG002.mat.cur.20211005.fasta
samtools faidx HG002.mat.cur.20211005.fasta SX | sed 's/SX/HG002#MAT#chrX/g' |\
  bgzip -@ 48 -c > hg002.chrX.fa.gz

~/tools/fastix/target/release/fastix-331c1159ea16625ee79d1a82522e800c99206834 -p 'HG002#PAT#' /lizardfs/erikg/T2T/liftover/split_chrY/chm13v2.0_chrY.fasta |\
  bgzip -@ 48 -c > hg002.chrY.fa.gz
```

Put HG002's chrX and chrY with the partitioned chrXs and chrYs.

```shell
# Prepare sequence order, with all references on the top
grep chr /lizardfs/erikg/HPRC/year1v2genbank/parts/chrS.pan.fa.fai | cut -f 1 > sequence_order.txt
grep chr /lizardfs/erikg/HPRC/year1v2genbank/parts/chrS.pan.fa.fai -v | cut -f 1 >> sequence_order.txt

cat \
  <(zcat hg002.chrX.fa.gz) \
  <(zcat hg002.chrY.fa.gz) \
  <(samtools faidx /lizardfs/erikg/HPRC/year1v2genbank/parts/chrS.pan.fa $(cat sequence_order.txt)) |\
   bgzip -@ 48 > chrS.pan+HG002chrXY.fa.gz
samtools faidx chrS.pan+HG002chrXY.fa.gz

rm sequence_order.txt
```

### Pangenome building

Apply `pggb` on the chromosome-partitioned HPRC dataset:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/graphs

num_of_haplotypes=$(cut -f 1,2 -d '#' /lizardfs/guarracino/chromosome_communities/assemblies/chrS.pan+HG002chrXY.fa.gz.fai | sort | uniq | wc -l)
sbatch -p highmem -c 48 --job-name sexpggb --wrap 'hostname; cd /scratch && /gnu/store/swnkjnc9wj6i1cl9iqa79chnf40r1327-pggb-0.2.0+640bf6b-5/bin/pggb -i /lizardfs/guarracino/chromosome_communities/assemblies/chrS.pan+HG002chrXY.fa.gz -o chrS.pan+HG002chrXY.s100k.l300k.p98.n'$num_of_haplotypes' -t 48 -s 100k -l 300k -p 98 -n '$num_of_haplotypes' -k 311 -G 13117,13219 -O 0.03 -T 48 -v -V chm13:#,grch38:#; mv /scratch/chrS.pan+HG002chrXY.s100k.l300k.p98.n'$num_of_haplotypes' /lizardfs/guarracino/chromosome_communities/graphs';
```

### Untangling

```shell
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrS.pan+HG002chrXY.s100k.l300k.p98.n93/chrS.pan+HG002chrXY.fa.gz.73e7992.4030258.b8e2fe5.smooth.fix.og
prefix=$(basename $path_input_og .og)

run_odgi=/home/guarracino/tools/odgi/bin/odgi-694948ccf31e7b565449cc056668e9dcc8cc0a3e
path_fasta_fai=/lizardfs/guarracino/chromosome_communities/assemblies/chrS.pan+HG002chrXY.fa.gz.fai

# All references and emit cut points
for refpattern in HG002; do
  path_targets_txt=/lizardfs/guarracino/chromosome_communities/untangle/$refpattern.target_paths.txt
  grep $refpattern $path_fasta_fai | cut -f 1 > $path_targets_txt
    
  for e in 5000 50000 100000; do
    for m in 500 1000 10000; do
      echo "-e $e -m $m"
    
      path_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$refpattern#SEX.e$e.m$m.bed.gz
      path_cut_points_txt=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$refpattern#SEX.e$e.m$m.cut_points.txt
    
      if [[ ! -s ${path_cut_points_txt} ]]; then
        sbatch -p workers -c 24 --job-name sexuntangle --wrap '\time -v '$run_odgi' untangle -t 24 -P -i '$path_input_og' -R '$path_targets_txt' -e '$e' -m '$m' --cut-points-output '$path_cut_points_txt' -j 0 -n 1 | pigz -c > '$path_bed_gz';'
      fi;
    done
  done
done

# Single reference by using the same cut points
for refpattern in HG002; do
  path_targets_txt=/lizardfs/guarracino/chromosome_communities/untangle/$refpattern.target_paths.txt
  
  for e in 5000 50000 100000; do
    for m in 500 1000 10000; do
      echo "-e $e -m $m"
        
      cat $path_targets_txt | while read ref; do
        echo $ref
            
        path_ref_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$ref.e$e.m$m.j0.n100.bed.gz
        if [[ ! -s ${path_ref_bed_gz} ]]; then
          path_cut_points_txt=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$refpattern#SEX.e$e.m$m.cut_points.txt
                    
          sbatch -p workers -c 24 --job-name sexuntangle --wrap '\time -v '$run_odgi' untangle -t 24 -P -i '$path_input_og' -r '$ref' -e '$e' -m '$m' --cut-points-input '$path_cut_points_txt' -j 0 -n 100 | pigz -c > '$path_ref_bed_gz';'
        fi;
      done
    done
  done
done

# Grounding
mkdir -p /lizardfs/guarracino/chromosome_communities/untangle/grounded

for refpattern in HG002; do
  path_targets_txt=/lizardfs/guarracino/chromosome_communities/untangle/$refpattern.target_paths.txt

  for e in 5000 50000 100000; do
    for m in 500 1000 10000; do
      for j in 0 0.8; do
        j_str=$(echo $j | sed 's/\.//g')
        (seq 1 5; seq 10 10 50) | while read n; do 
          echo "-e $e -m $m -j $j -n $n"
                
          cat $path_targets_txt | while read ref; do
            echo -e "\t"$ref

            path_grounded_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.tsv.gz
            if [[ ! -s ${path_grounded_tsv_gz} ]]; then
              path_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$refpattern#SEX.e$e.m$m.bed.gz
              path_ref_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$ref.e$e.m$m.j0.n100.bed.gz
              
              # Grounding
              ( echo query query.begin query.end target target.begin target.end jaccard strand self.coverage ref ref.begin ref.end | tr ' ' '\t'
                join \
                  <(zcat $path_bed_gz | awk '{ print $1"_"$2, $0 }' | tr ' ' '\t' | sort -k 1,1) \
                  <(zcat $path_ref_bed_gz | awk -v j=$j -v n=$n '{ if($7 >= j && $10 <= n) {print $1"_"$2, $0} }' | tr ' ' '\t' | sort -k 1,1) | \
                tr ' ' '\t' | grep -v '^#' | cut -f 2- | cut -f -9,14-16 ) | tr ' ' '\t' | pigz -c > x.tsv.gz
              
              # Contigs overlapping (or close at least 100kbps to) a PAR
              # Note that chrX PAR is from chm13, not HG002
              ref_chr=$(echo $ref | rev | cut -f 1 -d '#' | rev)
              bedtools intersect \
                -a <(zcat x.tsv.gz | awk -v OFS="\t" '{print $10,$11,$12,$1, "", "+"}' | grep query -v | bedtools sort) \
                -b <(grep $ref_chr /lizardfs/guarracino/chromosome_communities/data/chm13_hg002.PARs.bed |\
                  bedtools sort |\
                  bedtools slop -b 100000 -g /lizardfs/guarracino/chromosome_communities/data/chm13_hg002.PARs.sizes |\
                  awk -v OFS='\t' -v ref=$ref '{print(ref,$2,$3)}') | \
                #awk '$3-$2+1>=100000' | \
                cut -f 4 | \
                #Remove references to avoid grepping everything later (with zgrep -f)
                grep -v chr |\
                sort | uniq > $ref.tmp.txt
     
              # Add grounded.target column, re-add the references, and add annotation
              cat \
                <(zcat x.tsv.gz | head -n 1 | awk -v OFS='\t' '{print $0, "grounded.target"}') \
                <(zgrep x.tsv.gz -f $ref.tmp.txt | awk -v OFS='\t' -v ref=$ref '{print $0, ref}') \
                <(zcat x.tsv.gz | awk -v OFS='\t' -v ref=$ref '$1 ~ /chr/ { print $0, ref}' ) \
                <(grep $ref_chr /lizardfs/guarracino/chromosome_communities/data/chm13_hg002.PARs.approximate.bed | \
                  awk -v OFS='\t' -v ref=$ref '{print $1"#"$4,".",".",ref,".",".",".",".",".",".",$2,$3, ref}') |\
                pigz -c > $path_grounded_tsv_gz

              rm x.tsv.gz
            fi;

            ### Take only reliable blocks (flagged with "Hh" or "Hc", https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/README.md#components)???
          done
        done
      done
    done
  done
done
```

Plot:

```shell
for refpattern in HG002; do
  path_targets_txt=/lizardfs/guarracino/chromosome_communities/untangle/$refpattern.target_paths.txt

  for e in 5000 50000 100000; do
    for m in 500 1000 10000; do
      for j in 0 0.8; do
        j_str=$(echo $j | sed 's/\.//g')
        (seq 1 5; seq 10 10 50) | while read n; do 
          echo "-e $e -m $m -j $j -n $n"
    
          path_grounded_all_references_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$refpattern#chrSEX.e$e.m$m.j${j_str}.n$n.grounded.tsv.gz
          if [[ ! -s ${path_grounded_all_references_tsv_gz} ]]; then
            # Merge single reference results
            cat \
              <(zcat /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$refpattern*.e$e.m$m.j${j_str}.n$n.grounded.tsv.gz | head -n 1) \
              <(zcat /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$refpattern*.e$e.m$m.j${j_str}.n$n.grounded.tsv.gz | grep query -v) |\
              pigz -c > x.tsv.gz
            # Rename after to avoid getting itself with the previous '*' expansion
            mv x.tsv.gz $path_grounded_all_references_tsv_gz
          fi;

          Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_all.R $path_grounded_all_references_tsv_gz "-e $e -m $m -j $j -n $n" 155000000 800
        done
      done
    done
  done
done

mv /lizardfs/guarracino/chromosome_communities/untangle/grounded/*.pdf /lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/

# Merge
for refpattern in HG002; do
  path_targets_txt=/lizardfs/guarracino/chromosome_communities/untangle/$refpattern.target_paths.txt
  
  for e in 5000 50000 100000; do
    for m in 500 1000 10000; do
      for j in 0 0.8; do
        echo "-e $e -m $m -j $j"
            
        j_str=$(echo $j | sed 's/\.//g')
        PDFs=$((seq 1 5; seq 10 10 50) | while read n; do \
          echo /lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/$prefix.untangle.$refpattern#chrSEX.e$e.m$m.j${j_str}.n$n.grounded.tsv.gz.pdf
        done | tr '\n' ' ')
        #echo $PDFs

        /gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite $PDFs /lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/$prefix.untangle.$refpattern#chrSEX.e$e.m$m.j${j_str}.merged.grounded.tsv.gz.pdf
      done
    done
  done
done

#PAR1/2/3
# https://link.springer.com/article/10.1007/s10142-013-0323-6/figures/1
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

# sex-contigs
path_input_gfa=/lizardfs/guarracino/chromosome_communities/graphs/chrS.pan+HG002chrXY.s100k.l300k.p98.n93/chrS.pan+HG002chrXY.fa.gz.73e7992.4030258.b8e2fe5.smooth.fix.gfa
path_chm13_vcf_gz=/lizardfs/guarracino/chromosome_communities/graphs/chrS.pan+HG002chrXY.s100k.l300k.p98.n93/chrS.pan+HG002chrXY.fa.gz.73e7992.4030258.b8e2fe5.smooth.fix.chm13.vcf.gz
#path_grch38_vcf_gz=/lizardfs/guarracino/chromosome_communities/graphs/chrS.pan+HG002chrXY.s100k.l300k.p98.n93/chrS.pan+HG002chrXY.fa.gz.73e7992.4030258.b8e2fe5.smooth.fix.grch38.vcf.gz
path_hg002_vcf_gz=/lizardfs/guarracino/chromosome_communities/graphs/chrS.pan+HG002chrXY.s100k.l300k.p98.n93/chrS.pan+HG002chrXY.fa.gz.73e7992.4030258.b8e2fe5.smooth.fix.hg002.vcf.gz
sbatch -p workers -c 48 --job-name vgchm13 --wrap '\time -v vg deconstruct -P chm13 -H '?' -e -a -t 48 '$path_input_gfa' | bgzip -@ 48 -c > '$path_chm13_vcf_gz' && tabix '$path_chm13_vcf_gz
#sbatch -p workers -c 48 --job-name vggrch38 --wrap '\time -v vg deconstruct -P grch38 -H '?' -e -a -t 48 '$path_input_gfa' | bgzip -@ 48 -c > '$path_grch38_vcf_gz' && tabix '$path_grch38_vcf_gz
sbatch -p workers -c 48 --job-name vghg002 --wrap '\time -v vg deconstruct -P HG002 -H '?' -e -a -t 48 '$path_input_gfa' | bgzip -@ 48 -c > '$path_hg002_vcf_gz' && tabix '$path_hg002_vcf_gz


# In the VCF there are variants with all samples having missing genotype!
#num_samples=`bcftools query -l $PATH_VCF_GZ | wc -l`
#num_miss_gen=$(echo $num_samples - 1 | bc)
#--max-missing-count $num_miss_gen \
```
