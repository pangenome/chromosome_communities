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

(seq 13 15; seq 21 22) | while read f; do
  /gnu/store/d6ajy0ajxkbb22gkgi32g5c4jy8rs8bv-wfmash-0.6.0+948f168-21/bin/wfmash \
  /lizardfs/guarracino/chromosome_communities/assemblies/chm13.chr$f.fa.gz /lizardfs/erikg/HPRC/year1v2genbank/parts/chr$f.pan.fa -s 50k -l 150k -p 90 -n 1 -t 48 -m > chr$f.vs.chm13.chr$f.s50k.l150k.p90.n1.paf
done
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
