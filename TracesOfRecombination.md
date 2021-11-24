# Traces of recombination

```
(echo pggb wfmash seqwish smoothxg odgi gfaffix | tr ' ' '\n') | while read tool; do ls -l $(which $tool); done | cut -f 13 -d ' '
/gnu/store/2rwrch6gc5r5pazikhfc25j2am2rh22a-pggb-0.2.0+531f85f-1/bin/pggb
/gnu/store/nkfg1wg76zqaig43qgslkwcag9rb9fzz-wfmash-0.6.0+e9a5b02-17/bin/wfmash
/gnu/store/1v50lsaj5wjblznm3r8f83ccsw084wkr-seqwish-0.7.1+2ab95d7-6/bin/seqwish
/gnu/store/hyww54rrsga53xfr9gaixk08jq387zhj-smoothxg-0.6.0+d39280b-9/bin/smoothxg
/gnu/store/aw7x195pl3v83jk2ncni8xv1h82gxmwx-odgi-0.6.2+b04c7a9-11/bin/odgi
/gnu/store/31b6l8wxb6qgm4i8447233jv0z8bra5h-gfaffix-0.1.2.2/bin/gfaffix

vg
vg: variation graph tool, version v1.36.0 "Cibottola"
```

Create the main folder.

```
mkdir -p /lizardfs/guarracino/HPRC/chromosome_communities/
cd /lizardfs/guarracino/HPRC/chromosome_communities/
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
