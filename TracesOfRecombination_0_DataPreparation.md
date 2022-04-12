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
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/assemblies/year1_f1_assembly_v2_genbank/HG002.maternal.f1_assembly_v2_genbank.fa.gz
wget -c https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/assemblies/year1_f1_assembly_v2_genbank/HG002.paternal.f1_assembly_v2_genbank.fa.gz

for hap in maternal paternal; do
  echo $hap
  gunzip HG002.$hap.f1_assembly_v2_genbank.fa.gz -c | bgzip -@ 48 -c > HG002.$hap.fa.gz
  samtools faidx HG002.$hap.fa.gz
done
rm HG002.*.f1_assembly_v2_genbank.fa.gz

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
  reffa=/lizardfs/guarracino/chromosome_communities/assemblies/hg002-prox.renamed.$hap.fna.gz
  
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
ls partitioning/HG002.*.split.vs.hg002-prox.*.paf

for hap in mat pat; do
  echo $hap

  cat partitioning/HG002.$hap.split.vs.hg002-prox.$hap.paf | \
   awk '{ print $1,$11,$0 }' | tr ' ' '\t' |  sort -n -r -k 1,2 | \
   awk '$1 != last { print; last = $1; }' | awk -v OFS='\t' '{print($1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15)}' > partitioning/HG002.$hap.split.vs.hg002-prox.$hap.recovered.paf
done


## Subset by acrocentric chromosome
( seq 13 15; seq 21 22 ) | while read f; do  
  for hap in mat pat; do
    echo $f $hap
    hapfa=/lizardfs/guarracino/chromosome_communities/assemblies/HG002.${hap}ernal.fa.gz
  
    awk '$6 ~ "chr'$f'.prox$"' $( echo partitioning/HG002.$hap.vs.hg002-prox.$hap.paf; echo partitioning/HG002.$hap.split.vs.hg002-prox.$hap.recovered.paf ) | cut -f 1 > partitioning/hg002.$f.$hap.contigs.txt
    
    samtools faidx $hapfa $(cat partitioning/hg002.$f.$hap.contigs.txt) >> hg002.chr$f.hifi.fa
  done
  echo "bgzip..."
  
  bgzip -@ 48 hg002.chr$f.hifi.fa
  samtools faidx hg002.chr$f.hifi.fa.gz
done
```

Check the partitioning with CHM13:

```shell
## Map against the CHM13 assembly (and manually check the output)
for hap in mat pat; do
  reffa=/lizardfs/guarracino/chromosome_communities/assemblies/chm13.fa.gz
  hapfa=/lizardfs/guarracino/chromosome_communities/assemblies/HG002.${hap}ernal.fa.gz 
  
  sbatch -c 48 -p workers --job-name HG002 --wrap 'hostname; cd /scratch; \time -v /gnu/store/d6ajy0ajxkbb22gkgi32g5c4jy8rs8bv-wfmash-0.6.0+948f168-21/bin/wfmash -t 48 -m -N -s 50k -p 90 -n 1 '$reffa' '$hapfa' > /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/HG002.'$hap'.vs.CHM13.n1.paf'
done

for hap in mat pat; do
  reffa=/lizardfs/guarracino/chromosome_communities/assemblies/chm13.fa.gz
  hapfa=/lizardfs/guarracino/chromosome_communities/assemblies/HG002.${hap}ernal.fa.gz 
  
  sbatch -c 48 -p workers --job-name HG002 --wrap 'hostname; cd /scratch; \time -v /gnu/store/d6ajy0ajxkbb22gkgi32g5c4jy8rs8bv-wfmash-0.6.0+948f168-21/bin/wfmash -t 48 -m -s 50k -p 90 '$reffa' '$hapfa' > /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/HG002.'$hap'.vs.CHM13.split.paf'
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
