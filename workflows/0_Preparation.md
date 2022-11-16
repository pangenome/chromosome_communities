# Preparation

## Tools

```shell
mkdir -p ~/tools $$ cd ~/tools

git clone --recursive https://github.com/ekg/wfmash.git
cd wfmash
git checkout ad8aebae1be96847839778af534866bc9545adb9
cmake -H. -Bbuild && cmake --build build -- -j 48
mv build/bin/wfmash build/bin/wfmash-ad8aebae1be96847839778af534866bc9545adb9
cd ..

git clone --recursive https://github.com/ekg/seqwish.git
cd seqwish
git checkout 706ef7e2640e38a75ae7435fb57f7a6c8e3ada2c
cmake -H. -Bbuild && cmake --build build -- -j 48
mv bin/seqwish bin/seqwish-706ef7e2640e38a75ae7435fb57f7a6c8e3ada2c
cd ..

git clone --recursive https://github.com/pangenome/smoothxg.git
cd smoothxg
git checkout b3f4578f37000922cdd193c2d183951f48f6e612
cmake -H. -Bbuild && cmake --build build -- -j 48
mv bin/smoothxg bin/smoothxg-b3f4578f37000922cdd193c2d183951f48f6e612
cd ..

git clone --recursive https://github.com/pangenome/odgi.git
cd odgi
git checkout 454197fa29b772050c3135d5de47c816ce38e62c
cmake -H. -Bbuild && cmake --build build -- -j 48
mv bin/odgi bin/odgi-454197fa29b772050c3135d5de47c816ce38e62c

# For the bug-fixed odgi flip
git checkout 0b21b3525fc2ed9305c7df2386475a008a9337bd
cmake -H. -Bbuild && cmake --build build -- -j 48
mv bin/odgi bin/odgi-0b21b3525fc2ed9305c7df2386475a008a9337bd
cd ....

git clone --recursive https://github.com/pangenome/pggb.git
cd pggb
git checkout a4a6668d9ece42c80ce69dc354f0cb59a849286f
sed 's,"$fmt" wfmash,"$fmt" ~/tools/wfmash/build/bin/wfmash-ad8aebae1be96847839778af534866bc9545adb9,g' pggb -i
sed 's,"$fmt" seqwish,"$fmt" ~/tools/seqwish/bin/seqwish-706ef7e2640e38a75ae7435fb57f7a6c8e3ada2c,g' pggb -i
sed 's,"$fmt" smoothxg,"$fmt" ~/tools/smoothxg/bin/smoothxg-b3f4578f37000922cdd193c2d183951f48f6e612,g' pggb -i
sed 's,"$fmt" odgi,"$fmt" ~/tools/odgi/bin/odgi-454197fa29b772050c3135d5de47c816ce38e62c,g' pggb -i
mv pggb pggb-a4a6668d9ece42c80ce69dc354f0cb59a849286f
cd ..

git clone --recursive https://github.com/ekg/fastix.git
cd fastix
git checkout 331c1159ea16625ee79d1a82522e800c99206834
cargo build --release
mv target/release/fastix target/release/fastix-331c1159ea16625ee79d1a82522e800c99206834
cd ..

(echo gfaffix | tr ' ' '\n') | while read tool; do ls -l $(which $tool); done | cut -f 13 -d ' '
/gnu/store/1lfrwkhdpp95l8hh8grv75yhssn6in3r-gfaffix-0.1.3/bin/gfaffix

vg version | head -n 1
vg version v1.40.0 "Suardi"


gephi + "givecolortonodes" plugin (https://gephi.org/plugins/#/plugin/givecolortonodes)
```

## Data

Clone the repository:

```shell
cd /lizardfs/guarracino/
git clone --recursive https://github.com/pangenome/chromosome_communities.git
mkdir /lizardfs/guarracino/chromosome_communities/assemblies
cd /lizardfs/guarracino/chromosome_communities/assemblies
```

GRCh38's missing/masked (Ns) regions:

```shell
python3 /lizardfs/guarracino/chromosome_communities/scripts/generate_masked_ranges.py \
  /lizardfs/guarracino/chromosome_communities/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_auto_XY.fa.gz \
  > /lizardfs/guarracino/chromosome_communities/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_auto_XY.Ns.bed
awk -v OFS='\t' '{print($0,"grch38_Ns#000000")}' GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_auto_XY.Ns.bed \
  > /lizardfs/guarracino/chromosome_communities/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_auto_XY.Ns.black.bed
awk -v OFS='\t' '{print($0,"grch38_Ns#FFFFFF")}' GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_auto_XY.Ns.bed \
  > /lizardfs/guarracino/chromosome_communities/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_auto_XY.Ns.white.bed
```

### Whole HPRC dataset

Prepare a single FASTA with all HPRCy1v2genbank samples (94 haplotypes):

```shell
sbatch -p workers -c 48 --job-name collect --wrap 'cat /lizardfs/erikg/HPRC/year1v2genbank/assemblies/*v2_genbank*fa | bgzip -@ 48 -c > /lizardfs/guarracino/chromosome_communities/assemblies/HPRCy1v2genbank.fa.gz; samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/HPRCy1v2genbank.fa.gz'
```

### HG002

Prepare verkko's HG002 contigs (HiFi+ONT-based):

```shell
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

# HG002 has low quality, so no pq-contigs
mkdir -p partitioning

for hap in mat pat; do
  samtools faidx hg002-prox.renamed.fna $(grep $hap hg002-prox.renamed.fna.fai -i | cut -f 1) | bgzip -@ 48 -c > hg002-prox.renamed.$hap.fna.gz
  samtools faidx hg002-prox.renamed.$hap.fna.gz
done
```

Check manually the partitioning with CHM13:

```shell
RUN_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-ad8aebae1be96847839778af534866bc9545adb9

## Map against the CHM13 assembly (and manually check the output)
for hap in mat pat; do
  reffa=/lizardfs/guarracino/chromosome_communities/assemblies/chm13.fa.gz
  hapfa=/lizardfs/guarracino/chromosome_communities/assemblies/HG002.${hap}ernal.fa.gz 
  
  sbatch -c 48 -p workers --job-name HG002 --wrap "hostname; cd /scratch; \time -v $RUN_WFMASH -t 48 -m -N -s 50k -p 90 -n 1 $reffa $hapfa > /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/HG002.${hap}.vs.CHM13.n1.paf"
done

for hap in mat pat; do
  reffa=/lizardfs/guarracino/chromosome_communities/assemblies/chm13.fa.gz
  hapfa=/lizardfs/guarracino/chromosome_communities/assemblies/HG002.${hap}ernal.fa.gz 
  
  sbatch -c 48 -p workers --job-name HG002 --wrap "hostname; cd /scratch; \time -v $RUN_WFMASH -t 48 -m -s 50k -p 90 $reffa $hapfa > /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/HG002.${hap}.vs.CHM13.split.paf"
done
```

Prepare HG002-bakeoff contigs:

```shell
cd /lizardfs/guarracino/chromosome_communities/assemblies/

/home/guarracino/tools/fastix/target/release/fastix-331c1159ea16625ee79d1a82522e800c99206834 -p 'HG002-bakeoff#MAT#' <(zcat HG002.mat.cur.20211005.fasta.gz ) > hg002-bakeoff.mat.fa && samtools faidx hg002-bakeoff.mat.fa
/home/guarracino/tools/fastix/target/release/fastix-331c1159ea16625ee79d1a82522e800c99206834 -p 'HG002-bakeoff#PAT#' <(zcat HG002.pat.cur.20211005.fasta.gz ) > hg002-bakeoff.pat.fa && samtools faidx hg002-bakeoff.pat.fa
```

### HG01978 (verkko)

Prepare verkko's HG01978 contigs (HiFi+ONT-based):

```shell
### Get the "assembly.{m,p}at.scaff.fa.gz" files from Globus (t2t-share/HG01978/verkko-asm/v1 folder)

# Rename contigs
zcat assembly.mat.scaff.fa.gz | sed 's/^>mat/>HG01978#MAT#mat/g' > hg01978.mat.fa && samtools faidx hg01978.mat.fa
zcat assembly.pat.scaff.fa.gz | sed 's/^>pat/>HG01978#PAT#pat/g' > hg01978.pat.fa && samtools faidx hg01978.pat.fa
```

### Chromosome partitioning

Map contigs against the references:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/
cd /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/

REFERENCES_FASTA=/lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13+grch38.fa
RUN_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-ad8aebae1be96847839778af534866bc9545adb9

(ls /lizardfs/erikg/HPRC/year1v2genbank/assemblies/*v2_genbank*fa; \
  ls /lizardfs/guarracino/chromosome_communities/assemblies/hg01978.*at.fa; \
  ls /lizardfs/guarracino/chromosome_communities/assemblies/hg002-bakeoff.*at.fa) | while read FASTA; do
  HAPLOTYPE=$(basename $FASTA .fa | cut -f 1,2 -d '.');
  echo $HAPLOTYPE
  
  PAF=/lizardfs/guarracino/chromosome_communities/assemblies/partitioning/$HAPLOTYPE.vs.ref.paf
  sbatch -p workers -c 12 --wrap "$RUN_WFMASH -t 12 -m -N -s 50k -l 150k -p 90 -H 0.001 $REFERENCES_FASTA $FASTA > $PAF"
done
```

Collect unmapped contigs and remap them in split mode:

```shell
(ls /lizardfs/erikg/HPRC/year1v2genbank/assemblies/*v2_genbank*fa; \
  ls /lizardfs/guarracino/chromosome_communities/assemblies/hg01978.*at.fa; \
  ls /lizardfs/guarracino/chromosome_communities/assemblies/hg002-bakeoff.*at.fa) | while read FASTA; do
  HAPLOTYPE=$(basename $FASTA .fa | cut -f 1,2 -d '.');
  echo $HAPLOTYPE
  
  UNALIGNED=/lizardfs/guarracino/chromosome_communities/assemblies/partitioning/$HAPLOTYPE.unaligned
  
  PAF=/lizardfs/guarracino/chromosome_communities/assemblies/partitioning/$HAPLOTYPE.vs.ref.paf
  comm -23 <(cut -f 1 $FASTA.fai | sort) <(cut -f 1 $PAF | sort) > $UNALIGNED.txt
  if [[ $(wc -l $UNALIGNED.txt | cut -f 1 -d\ ) != 0 ]];
  then 
    samtools faidx $FASTA $(tr '\n' ' ' < $UNALIGNED.txt) > $UNALIGNED.fa
    samtools faidx $UNALIGNED.fa
    sbatch -p workers -c 12 --wrap "$RUN_WFMASH -t 12 -m -s 50k -l 150k -p 90 -H 0.001 $REFERENCES_FASTA $UNALIGNED.fa > $UNALIGNED.split.vs.ref.paf"
  fi
done
```

Collect our best mapping for each of our attempted split rescues:

```shell
ls *.unaligned.split.vs.ref.paf | while read PAF; do
  cat $PAF | awk -v OFS='\t' '{ print $1,$11,$0 }' | sort -n -r -k 1,2 | \
    awk -v OFS='\t' '$1 != last { print($3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15); last = $1; }'
done > rescues.paf
```

More precise mapping information (pq-arms):

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/pq_info/
cd /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/pq_info/

# Prepare p/q-arms coordinates
sed 's/chr/chm13#chr/g' /lizardfs/guarracino/chromosome_communities/data/chm13.centromeres.approximate.bed | \
  bedtools sort | \
  bedtools complement \
    -i - \
    -g <(cut -f 1,2 /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa.fai | grep 'chrM' -v | sort) \
  > tmp.bed
sed 's/chr/grch38#chr/g' /lizardfs/guarracino/chromosome_communities/data/grch38.centromere.bed | \
  bedtools sort | \
  bedtools complement \
    -i - \
    -g <(cut -f 1,2 /lizardfs/erikg/HPRC/year1v2genbank/assemblies/grch38.fa.fai| grep 'chrM\|chrEBV' -v | grep '_' -v | sort) \
  >> tmp.bed
# Take odd rows
(echo -e "#chrom\tstart\tend"; sed -n 1~2p tmp.bed) > p_arms.bed
# Take even rows
(echo -e "#chrom\tstart\tend"; sed -n 2~2p tmp.bed) > q_arms.bed
rm tmp.bed

# Classify partitioned contigs with pq information
ls /lizardfs/erikg/HPRC/year1v2genbank/approx_mappings/*.vs.ref.paf | while read PAF; do
  HAPLOTYPE=$(basename $PAF .paf | cut -f 1,2 -d '.');
  echo $HAPLOTYPE
  
  # This instruction helps with PAF files made by enabling the mapping split in wfmash (no `-N` option).
  # With unsplit PAF files, this instruction has no effect
  cat $PAF | awk -v OFS='\t' '{ print $1,$11,$0 }' | sort -n -r -k 1,2 | \
    awk -v OFS='\t' '$1 != last { print($3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15); last = $1; }' | \
    awk -v OFS='\t' '{print($6,$8,$9,$1)}' > $HAPLOTYPE.vs.ref.bed

  bedtools intersect \
    -a <(cat $HAPLOTYPE.vs.ref.bed | bedtools sort) \
    -b <(cat p_arms.bed | bedtools sort) | cut -f 4 | sort > $HAPLOTYPE.vs.ref.p_contigs.txt
  bedtools intersect \
    -a <(cat $HAPLOTYPE.vs.ref.bed | bedtools sort) \
    -b <(cat q_arms.bed | bedtools sort) | cut -f 4 | sort > $HAPLOTYPE.vs.ref.q_contigs.txt
    
  grep -f <(comm -23 $HAPLOTYPE.vs.ref.p_contigs.txt $HAPLOTYPE.vs.ref.q_contigs.txt) $PAF | \
    awk -v OFS='\t' '{print($1,$6"_p")}' > $HAPLOTYPE.partitioning_with_pq.tsv
  grep -f <(comm -13 $HAPLOTYPE.vs.ref.p_contigs.txt $HAPLOTYPE.vs.ref.q_contigs.txt) $PAF | \
    awk -v OFS='\t' '{print($1,$6"_q")}' >> $HAPLOTYPE.partitioning_with_pq.tsv
  grep -f <(comm -12 $HAPLOTYPE.vs.ref.p_contigs.txt $HAPLOTYPE.vs.ref.q_contigs.txt) $PAF | \
    awk -v OFS='\t' '{print($1,$6"_pq")}' >> $HAPLOTYPE.partitioning_with_pq.tsv
    
  rm $HAPLOTYPE.vs.ref.bed $HAPLOTYPE.vs.ref.p_contigs.txt $HAPLOTYPE.vs.ref.q_contigs.txt
done
```

### Chromosome partitioning (alternative approach)

Put the CHM13 (v2.0) and GRCh38 references together:

```shell
cd /lizardfs/guarracino/chromosome_communities/assemblies/

cat \
  <( /home/guarracino/tools/fastix/target/release/fastix-331c1159ea16625ee79d1a82522e800c99206834 -p 'chm13#' <(zcat /lizardfs/erikg/human/chm13v2.0.fa.gz) ) \
  <( zcat GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_auto_XY.fa.gz) \
  <( samtools faidx /lizardfs/guarracino/pggb-paper/hsapiens/assemblies/grch38.fa.gz $(echo "grch38#chrM" "grch38#chrEBV") ) | \
  bgzip -c -@ 48 > chm13v2+grch38masked.fa.gz
  samtools faidx chm13v2+grch38masked.fa.gz
```

Map contigs against the references:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/
cd /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/

REFERENCES_FASTA=/lizardfs/guarracino/chromosome_communities/assemblies/chm13v2+grch38masked.fa.gz
RUN_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-ad8aebae1be96847839778af534866bc9545adb9

(ls /lizardfs/erikg/HPRC/year1v2genbank/assemblies/*v2_genbank*fa; \
  ls /lizardfs/guarracino/chromosome_communities/assemblies/hg01978.*at.fa; \
  ls /lizardfs/guarracino/chromosome_communities/assemblies/hg002-bakeoff.*at.fa) | while read FASTA; do
  HAPLOTYPE=$(basename $FASTA .fa | cut -f 1,2 -d '.');
  echo $HAPLOTYPE
  
  PAF=/lizardfs/guarracino/chromosome_communities/assemblies/partitioning/$HAPLOTYPE.vs.refs.paf
  sbatch -p workers -c 12 --wrap "$RUN_WFMASH -t 12 -m -s 5k -p 90 -H 0.001 $REFERENCES_FASTA $FASTA > $PAF"
done
```

Get the best target for each contig:


```shell
cd /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/

(ls /lizardfs/erikg/HPRC/year1v2genbank/assemblies/*v2_genbank*fa; \
  ls /lizardfs/guarracino/chromosome_communities/assemblies/hg01978.*at.fa; \
  ls /lizardfs/guarracino/chromosome_communities/assemblies/hg002-bakeoff.*at.fa) | while read FASTA; do
  HAPLOTYPE=$(basename $FASTA .fa | cut -f 1,2 -d '.');
  echo $HAPLOTYPE
  
  PAF=/lizardfs/guarracino/chromosome_communities/assemblies/partitioning/$HAPLOTYPE.vs.refs.paf

  # We remove the reference prefixes because we are interested in the chromosome
  python3 /lizardfs/guarracino/chromosome_communities/scripts/partition_by_weighted_mappings.py <(sed 's/chm13#//g;s/grch38#//g' $PAF) > $HAPLOTYPE.vs.refs.partitions.tsv
done
```

Get acrocentric contigs:

```shell
cd /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/

REFERENCES_FASTA=/lizardfs/guarracino/chromosome_communities/assemblies/chm13v2+grch38masked.fa.gz

rm chrACRO+refs.fa.gz
rm chrACRO.fa.gz
(seq 13 15; seq 21 22) | while read i; do
  echo chr$i
  
  samtools faidx $REFERENCES_FASTA $(grep chr$i $REFERENCES_FASTA.fai | cut -f 1) >> chrACRO+refs.fa
  
  zcat /lizardfs/guarracino/chromosome_communities/assemblies/hg002.chr$i.prox.fa.gz >> chrACRO+refs.fa
  zcat /lizardfs/guarracino/chromosome_communities/assemblies/hg002.chr$i.prox.fa.gz >> chrACRO.fa
  
  ls /lizardfs/erikg/HPRC/year1v2genbank/assemblies/*v2_genbank*fa | while read FASTA; do
    HAPLOTYPE=$(basename $FASTA .fa | cut -f 1,2 -d '.');
    echo $HAPLOTYPE
      
    samtools faidx $FASTA $(grep chr$i $HAPLOTYPE.vs.refs.partitions.tsv | cut -f 1) >> chrACRO+refs.fa
    samtools faidx $FASTA $(grep chr$i $HAPLOTYPE.vs.refs.partitions.tsv | cut -f 1) >> chrACRO.fa
  done
done
bgzip -@ 48 chrACRO+refs.fa && samtools faidx chrACRO+refs.fa.gz
bgzip -@ 48 chrACRO.fa && samtools faidx chrACRO.fa.gz



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
