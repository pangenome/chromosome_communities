# Untangling

## Acrocentric chromosomes

### Collect contigs running from the p-arm to the q-arm of the acrocentric chromosomes

Map contigs against the reference:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/assemblies/partitioning_pq/
cd /lizardfs/guarracino/chromosome_communities/assemblies/partitioning_pq/


REFERENCE_FASTA=/lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa
RUN_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-ad8aebae1be96847839778af534866bc9545adb9

(ls /lizardfs/erikg/HPRC/year1v2genbank/assemblies/*v2_genbank*fa; \
  ls /lizardfs/guarracino/chromosome_communities/assemblies/hg01978.*at.fa; \
  ls /lizardfs/guarracino/chromosome_communities/assemblies/hg002-bakeoff.*at.fa) | while read FASTA; do
  HAPLOTYPE=$(basename $FASTA .fa | cut -f 1,2 -d '.');
  echo $HAPLOTYPE
  
  PAF=/lizardfs/guarracino/chromosome_communities/assemblies/partitioning_pq/$HAPLOTYPE.vs.ref.paf
  sbatch -p workers -c 12 --wrap "$RUN_WFMASH -t 12 -m -s 50k -l 150k -p 90 -H 0.001 $REFERENCE_FASTA $FASTA > $PAF"
done

# HiFi HG002 contigs are partitioned with `wfmash -N` (it is more strict) because they are not filtered to be pq-contigs
ls /lizardfs/erikg/HPRC/year1v2genbank/assemblies/HG002*v2_genbank*fa | while read FASTA; do
  HAPLOTYPE=$(basename $FASTA .fa | cut -f 1,2 -d '.');
  echo $HAPLOTYPE
  
  PAF=/lizardfs/guarracino/chromosome_communities/assemblies/partitioning_pq/$HAPLOTYPE.vs.ref.p90.N.paf
  sbatch -p workers -c 12 --wrap "$RUN_WFMASH -t 12 -m -s 50k -l 150k -p 90 -N -H 0.001 $REFERENCE_FASTA $FASTA > $PAF"
done
```


Prepare the BED files for p-arms and q-arms (-/+ 1 Mbps):

```shell
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
sed -n 1~2p tmp.bed > p_arms.bed
# Take even rows
sed -n 2~2p tmp.bed > q_arms.bed

rm tmp.bed
```


Find the contigs which have mappings at least 1kbp-long in both the p-arm and the q-arm of the same chromosome:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/pq_contigs
cd /lizardfs/guarracino/chromosome_communities/pq_contigs

(seq 13 15; seq 21 22) | while read f; do
  ref=chm13#chr$f
  
  if [[ ! -s $ref.pq_contigs.1kbps.txt ]]; then
    echo $ref
    
    # "p-touching contigs" intersected "q-touching contigs"
    comm -12 \
      <(bedtools intersect \
        -a <(cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning_pq/*.paf | grep $ref | awk -v OFS='\t' '{print $6, $8, $9, $1, "", "+"}' | bedtools sort) \
        -b <(grep $ref /lizardfs/guarracino/chromosome_communities/assemblies/partitioning_pq/p_arms.bed | bedtools sort) | \
        awk '$3-$2+1>=1000' | \
        cut -f 4 | sort | uniq) \
      <(bedtools intersect \
        -a <(cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning_pq/*.paf | grep $ref | awk -v OFS='\t' '{print $6, $8, $9, $1, "", "+"}' | bedtools sort) \
        -b <(grep $ref /lizardfs/guarracino/chromosome_communities/assemblies/partitioning_pq/q_arms.bed | bedtools sort) | \
        awk '$3-$2+1>=1000' | \
        cut -f 4 | sort | uniq) > $ref.pq_contigs.1kbps.txt
  fi;
done

# Num. of contigs
#ls *txt | while read f; do echo $f; cat $f | wc -l; done
```


Prepare the FASTA files with the pq-contigs:

```shell
cd /lizardfs/guarracino/chromosome_communities/pq_contigs

PATH_HPRCY1_FA_GZ=/lizardfs/guarracino/chromosome_communities/assemblies/HPRCy1v2genbank.fa.gz

(seq 13 15; seq 21 22) | while read f; do
  ref=chm13#chr$f
  
  PATH_PQ_CONTIGS_FA_GZ=$ref.pq_contigs.1kbps.hg002prox.fa.gz
  if [[ ! -s $PATH_PQ_CONTIGS_FA_GZ ]]; then
    echo $ref
    
    cat \
      <(samtools faidx $PATH_HPRCY1_FA_GZ $(cat $ref.pq_contigs.1kbps.txt)) \
      <(zcat /lizardfs/guarracino/chromosome_communities/assemblies/hg002.chr$f.prox.fa.gz) | \
      bgzip -@ 48 -c > $PATH_PQ_CONTIGS_FA_GZ
    samtools faidx $PATH_PQ_CONTIGS_FA_GZ
  fi;
done

# Num. of contigs
#(seq 13 15; seq 21 22) | while read f; do echo -n "$f -> "; ref=chm13#chr$f; cat $ref.pq_contigs.1kbps.hg002prox.fa.gz.fai | wc -l; done
```


Put all together with/without both the human references:

```shell
cd /lizardfs/guarracino/chromosome_communities/pq_contigs

PATH_CHM13_FA=/lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa
PATH_GRCH38_FA_GZ=/lizardfs/guarracino/chromosome_communities/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_auto_XY.fa.gz

cat \
  <( samtools faidx $PATH_CHM13_FA $(grep 'chr13\|chr14\|chr15\|chr21\|chr22' $PATH_CHM13_FA.fai | cut -f 1) ) \
  <( samtools faidx $PATH_GRCH38_FA_GZ $(grep 'chr13\|chr14\|chr15\|chr21\|chr22' $PATH_GRCH38_FA_GZ.fai | grep '_' -v | cut -f 1) ) \
  <( zcat chm13#chr*.pq_contigs.1kbps.hg002prox.fa.gz ) | \
  bgzip -@ 48 -c > chrACRO+refs.pq_contigs.1kbps.hg002prox.fa.gz
samtools faidx chrACRO+refs.pq_contigs.1kbps.hg002prox.fa.gz
```


Put all together with both the human references plus all partitioned HG002 HiFi contigs:

```shell
cd /lizardfs/guarracino/chromosome_communities/pq_contigs

# Include long HG002 contigs that were not already taken and that map on the right of the rDNA array
bedtools complement \
  -i <( zgrep rDNA /lizardfs/guarracino/chromosome_communities/data/chm13.CentromeresAndTelomeres.CenSatAnnotation.txt.gz | cut -f 1,2,3 |  sed 's/chr/chm13#chr/g' | sort ) \
  -g <( cut -f 1,2 /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa.fai | grep 'chr13\|chr14\|chr15\|chr21\|chr22' | sort ) | \
  sed -n 2~2p > acrocentrics.prox.bed

comm -23 \
  <( bedtools intersect \
       -a <(cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning_pq/HG002.*.vs.ref.p90.N.paf | awk '$2 >= 300000' | awk -v OFS='\t' '{print $6, $8, $9, $1, "", "+"}' | bedtools sort) \
       -b acrocentrics.prox.bed | cut -f 4 | sort | uniq ) \
  <( cat chm13#chr*.pq_contigs.1kbps.hg002prox.fa.gz.fai | grep HG002 | cut -f 1 | sort | uniq ) | sort | uniq > HG002.contigs.txt

cat \
  <( samtools faidx $PATH_CHM13_FA $(grep 'chr13\|chr14\|chr15\|chr21\|chr22' $PATH_CHM13_FA.fai | cut -f 1) ) \
  <( samtools faidx $PATH_GRCH38_FA_GZ $(grep 'chr13\|chr14\|chr15\|chr21\|chr22' $PATH_GRCH38_FA_GZ.fai | grep '_' -v | cut -f 1) ) \
  <( zcat chm13#chr*.pq_contigs.1kbps.hg002prox.fa.gz ) \
  <( samtools faidx /lizardfs/erikg/HPRC/year1v2genbank/assemblies/HG002.paternal.f1_assembly_v2_genbank.fa $( grep '#1#' HG002.contigs.txt ) ) \
  <( samtools faidx /lizardfs/erikg/HPRC/year1v2genbank/assemblies/HG002.maternal.f1_assembly_v2_genbank.fa $( grep '#2#' HG002.contigs.txt ) ) | \
  bgzip -@ 48 -c > chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz
samtools faidx chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz

# Without references
cat \
  <( zcat chm13#chr*.pq_contigs.1kbps.hg002prox.fa.gz ) \
  <( samtools faidx /lizardfs/erikg/HPRC/year1v2genbank/assemblies/HG002.paternal.f1_assembly_v2_genbank.fa $( grep '#1#' HG002.contigs.txt ) ) \
  <( samtools faidx /lizardfs/erikg/HPRC/year1v2genbank/assemblies/HG002.maternal.f1_assembly_v2_genbank.fa $( grep '#2#' HG002.contigs.txt ) ) | \
  bgzip -@ 48 -c > chrACRO.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz
samtools faidx chrACRO.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz
```


### Collect contigs running from the p-arm to the q-arm of the acrocentric chromosomes (alternative approach)

Map contigs against the CHM13v2 reference:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/assemblies/partitioning_pq/s5k
cd /lizardfs/guarracino/chromosome_communities/assemblies/partitioning_pq/s5k

REFERENCE_FASTA=/lizardfs/erikg/human/chm13v2.0.fa.gz
RUN_WFMASH=/home/guarracino/tools/wfmash/build/bin/wfmash-ad8aebae1be96847839778af534866bc9545adb9

(ls /lizardfs/erikg/HPRC/year1v2genbank/assemblies/*v2_genbank*fa; \
  ls /lizardfs/guarracino/chromosome_communities/assemblies/hg01978.*at.fa; \
  ls /lizardfs/guarracino/chromosome_communities/assemblies/hg002-bakeoff.*at.fa) | while read FASTA; do
  HAPLOTYPE=$(basename $FASTA .fa | cut -f 1,2 -d '.');
  echo $HAPLOTYPE
  
  PAF=/lizardfs/guarracino/chromosome_communities/assemblies/partitioning_pq/s5k/$HAPLOTYPE.vs.ref.paf
  sbatch -p workers -c 12 --wrap "$RUN_WFMASH -t 12 -m -s 5k -p 90 -H 0.001 $REFERENCE_FASTA $FASTA > $PAF"
done
```


Prepare the BED files for p-arms and q-arms (-/+ 0.5 Mbps):

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/pq_contigs/s5k
cd /lizardfs/guarracino/chromosome_communities/pq_contigs/s5k

bedtools slop \
    -i <(sed 's/chr/chm13#chr/g' /lizardfs/guarracino/chromosome_communities/data/chm13.centromeres.approximate.bed | bedtools sort) \
    -g <(cut -f 1,2 /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa.fai | sort) \
    -b 500000 | \
  bedtools sort | \
  bedtools complement \
    -i - \
    -g <(cut -f 1,2 /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa.fai | sort) | \
  grep 'chr13\|chr14\|chr15\|chr21\|chr22' > tmp.bed
  
# Take odd rows
sed -n 1~2p tmp.bed > p_arms.bed
# Take even rows
sed -n 2~2p tmp.bed > q_arms.bed

rm tmp.bed
```


Find the contigs which have mappings at least 1kbp-long in both the p-arm and the q-arm of the same chromosome:

```shell
cd /lizardfs/guarracino/chromosome_communities/pq_contigs/s5k

# "p-touching contigs" intersected "q-touching contigs"
comm -12 \
  <(bedtools intersect \
    -a <(cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning_pq/s5k/*.paf | awk -v OFS='\t' '{print $6, $8, $9, $1, "", "+"}' | sed 's/^/chm13#/g' | bedtools sort) \
    -b <(cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning_pq/p_arms.bed | bedtools sort) | \
    awk '$3-$2+1>=1000' | \
    cut -f 4 | sort | uniq) \
  <(bedtools intersect \
    -a <(cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning_pq/s5k/*.paf | awk -v OFS='\t' '{print $6, $8, $9, $1, "", "+"}' | sed 's/^/chm13#/g' | bedtools sort) \
    -b <(cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning_pq/q_arms.bed | bedtools sort) | \
    awk '$3-$2+1>=1000' | \
    cut -f 4 | sort | uniq) > all.pq_contigs.1kbps.txt


# Partition pq-contigs by chromosome
(seq 13 15; seq 21 22) | while read f; do
  ref=chm13#chr$f
  
  if [[ ! -s $ref.pq_contigs.1kbps.txt ]]; then
    echo $ref

    comm -12 \
      <(sort all.pq_contigs.1kbps.txt) \
      <(cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/*.partitions.tsv | awk -v target=chr$f -v OFS='\t' '$2 == target {print $1}' | sort) \
      > $ref.pq_contigs.1kbps.txt
  fi;
done

# Num. of contigs
ls *txt | while read f; do echo $f; cat $f | wc -l; done
```


Prepare the FASTA files with the pq-contigs:

```shell
cd /lizardfs/guarracino/chromosome_communities/pq_contigs/s5k

PATH_HPRCY1_FA_GZ=/lizardfs/guarracino/chromosome_communities/assemblies/HPRCy1v2genbank.fa.gz

(seq 13 15; seq 21 22) | while read f; do
  ref=chm13#chr$f
  
  PATH_PQ_CONTIGS_FA_GZ=$ref.pq_contigs.1kbps.hg002prox.fa.gz
  if [[ ! -s $PATH_PQ_CONTIGS_FA_GZ ]]; then
    echo $ref
    
    # HiFi-based assemblies + hg01978 verkko + hg002-bakeoff + hg002 verkko
    cat \
      <(samtools faidx $PATH_HPRCY1_FA_GZ $(comm -12 <(cut -f 1 "$PATH_HPRCY1_FA_GZ".fai | sort) <(sort $ref.pq_contigs.1kbps.txt))) \
      <(samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/hg01978.mat.fa $(comm -12 <(cut -f 1 /lizardfs/guarracino/chromosome_communities/assemblies/hg01978.mat.fa.fai | sort) <(sort $ref.pq_contigs.1kbps.txt))) \
      <(samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/hg01978.pat.fa $(comm -12 <(cut -f 1 /lizardfs/guarracino/chromosome_communities/assemblies/hg01978.pat.fa.fai | sort) <(sort $ref.pq_contigs.1kbps.txt))) \
      <(samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/hg002-bakeoff.mat.fa $(comm -12 <(cut -f 1 /lizardfs/guarracino/chromosome_communities/assemblies/hg002-bakeoff.mat.fa.fai | sort) <(sort $ref.pq_contigs.1kbps.txt))) \
      <(samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/hg002-bakeoff.pat.fa $(comm -12 <(cut -f 1 /lizardfs/guarracino/chromosome_communities/assemblies/hg002-bakeoff.pat.fa.fai | sort) <(sort $ref.pq_contigs.1kbps.txt))) \
      <(zcat /lizardfs/guarracino/chromosome_communities/assemblies/hg002.chr$f.prox.fa.gz) | \
      bgzip -@ 48 -c > $PATH_PQ_CONTIGS_FA_GZ
    samtools faidx $PATH_PQ_CONTIGS_FA_GZ
  fi;
done

# Num. of contigs
(seq 13 15; seq 21 22) | while read f; do echo -n "$f -> "; ref=chm13#chr$f; cat $ref.pq_contigs.1kbps.hg002prox.fa.gz.fai | wc -l; done
```


Put all together with both the human references:

```shell
cd /lizardfs/guarracino/chromosome_communities/pq_contigs/s5k

PATH_CHM13_FA_GZ=/lizardfs/erikg/human/chm13v2.0.fa.gz
PATH_GRCH38_FA_GZ=/lizardfs/guarracino/chromosome_communities/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_auto_XY.fa.gz

cat \
  <( samtools faidx $PATH_CHM13_FA_GZ $(grep 'chr13\|chr14\|chr15\|chr21\|chr22' $PATH_CHM13_FA_GZ.fai | cut -f 1) | sed 's/^>/>chm13#/g' ) \
  <( samtools faidx $PATH_GRCH38_FA_GZ $(grep 'chr13\|chr14\|chr15\|chr21\|chr22' $PATH_GRCH38_FA_GZ.fai | grep '_' -v | cut -f 1) ) \
  <( zcat chm13#chr*.pq_contigs.1kbps.hg002prox.fa.gz ) | \
  bgzip -@ 48 -c > chrACRO+refs.pq_contigs.1kbps.hg002prox.fa.gz
samtools faidx chrACRO+refs.pq_contigs.1kbps.hg002prox.fa.gz
```


Put all together with both the human references plus all partitioned HG002 HiFi contigs:

```shell
cd /lizardfs/guarracino/chromosome_communities/pq_contigs/s5k

# Include long HG002 contigs that were not already taken and that map on the right of the rDNA array
bedtools complement \
  -i <( zgrep rDNA /lizardfs/guarracino/chromosome_communities/data/chm13.CentromeresAndTelomeres.CenSatAnnotation.txt.gz | cut -f 1,2,3 |  sed 's/chr/chm13#chr/g' | sort ) \
  -g <( cut -f 1,2 /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa.fai | grep 'chr13\|chr14\|chr15\|chr21\|chr22' | sort ) | \
  sed -n 2~2p > acrocentrics.prox.bed

comm -12 \
  <( bedtools intersect \
       -a <(cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning_pq/s5k/HG002.*.vs.ref.paf | awk '$2 >= 300000' | awk -v OFS='\t' '{print "chm13#"$6, $8, $9, $1, "", "+"}' | bedtools sort) \
       -b acrocentrics.prox.bed | cut -f 4 | sort | uniq ) \
  <( cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/HG002.*.partitions.tsv | grep 'chr13\|chr14\|chr15\|chr21\|chr22' | sort | cut -f 1 ) | sort | uniq > HG002.tmp.contigs.txt

# To avoid duplicates
comm -23 HG002.tmp.contigs.txt <(cat chm13#chr*.pq_contigs.1kbps.hg002prox.fa.gz.fai | grep HG002 | cut -f 1 | sort | uniq) > HG002.contigs.txt
rm HG002.tmp.contigs.txt

cat \
  <( samtools faidx $PATH_CHM13_FA_GZ $(grep 'chr13\|chr14\|chr15\|chr21\|chr22' $PATH_CHM13_FA_GZ.fai | cut -f 1) | sed 's/^>/>chm13#/g' ) \
  <( samtools faidx $PATH_GRCH38_FA_GZ $(grep 'chr13\|chr14\|chr15\|chr21\|chr22' $PATH_GRCH38_FA_GZ.fai | grep '_' -v | cut -f 1) ) \
  <( zcat chm13#chr*.pq_contigs.1kbps.hg002prox.fa.gz ) \
  <( samtools faidx /lizardfs/erikg/HPRC/year1v2genbank/assemblies/HG002.paternal.f1_assembly_v2_genbank.fa $( grep '#1#' HG002.contigs.txt ) ) \
  <( samtools faidx /lizardfs/erikg/HPRC/year1v2genbank/assemblies/HG002.maternal.f1_assembly_v2_genbank.fa $( grep '#2#' HG002.contigs.txt ) ) | \
  bgzip -@ 48 -c > chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz
samtools faidx chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz
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

RUN_PGGB=/home/guarracino/tools/pggb/pggb-a4a6668d9ece42c80ce69dc354f0cb59a849286f

num_of_haplotypes=$(cat /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.fa.gz.fai | wc -l)
num_of_haplotypes_plus_a_bit=$(echo "$num_of_haplotypes + 10" | bc)  # 5 mat acros + 5 pat acros
sbatch -p workers -c 48 --job-name acropggb --wrap "hostname; cd /scratch && $RUN_PGGB -i /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz -o chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n${num_of_haplotypes_plus_a_bit} -t 48 -s 50k -l 250k -p 98 -F 0.001 -n ${num_of_haplotypes_plus_a_bit} -k 311 -G 13117,13219 -O 0.03 -T 48 -V chm13:#,grch38:#; mv /scratch/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n${num_of_haplotypes_plus_a_bit} /lizardfs/guarracino/chromosome_communities/graphs";

# GFA with only chm13's paths (for faster loading with gfaestus)
f=chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final
grep '^H\|^S\|^L' $f.gfa > $f.only_chm13_paths.gfa
grep 'chm13#chr' $f.gfa >> $f.only_chm13_paths.gfa

# GFA with only chm13 and grch38's paths (for faster loading with gfaestus)
grep '^H\|^S\|^L' $f.gfa > $f.only_chm13+grch38_paths.gfa
grep 'chm13#chr\|grch38#chr' $f.gfa >> $f.only_chm13+grch38_paths.gfa

# -x 200 -G 20 -I 10000 -l 1000
sbatch -c 48 -p workers --wrap "cd /scratch && /home/guarracino/tools/odgi/bin/odgi-454197fa29b772050c3135d5de47c816ce38e62c layout -i /lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.og -o /lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.og.x200.G20.I10000.l1000.lay -T /lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.og.x200.G20.I10000.l1000.lay.tsv -x 200 -G 20 -I 10000 -l 1000 -t 48 -P"

#todo to finish
seqtk cutN -n 10000 -g /lizardfs/guarracino/chromosome_communities/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_auto_XY.fa.gz > /lizardfs/guarracino/chromosome_communities/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_auto_XY.cutN.n10000.bed
seqtk cutN -n 1000  -g /lizardfs/guarracino/chromosome_communities/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_auto_XY.fa.gz > /lizardfs/guarracino/chromosome_communities/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_auto_XY.cutN.n1000.bed
seqtk cutN -n 100   -g /lizardfs/guarracino/chromosome_communities/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_auto_XY.fa.gz > /lizardfs/guarracino/chromosome_communities/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_auto_XY.cutN.n100.bed
seqtk cutN -n 10    -g /lizardfs/guarracino/chromosome_communities/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_auto_XY.fa.gz > /lizardfs/guarracino/chromosome_communities/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_auto_XY.cutN.n10.bed
seqtk cutN -n 1     -g /lizardfs/guarracino/chromosome_communities/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_auto_XY.fa.gz > /lizardfs/guarracino/chromosome_communities/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_auto_XY.cutN.n1.bed

# Reference-free pq-contigs acrocentric graph
#num_of_haplotypes=$(cat /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.fa.gz.fai | wc -l)
num_of_haplotypes=10
sbatch -p workers -c 48 --job-name rfreeacropggb --wrap "hostname; cd /scratch && $RUN_PGGB -i /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz -o chrACRO.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n${num_of_haplotypes_plus_a_bit} -t 48 -s 50k -l 250k -p 98 -F 0.001 -n ${num_of_haplotypes_plus_a_bit} -k 311 -G 13117,13219 -O 0.03 -T 48 -V chm13:#,grch38:#; mv /scratch/chrACRO.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n${num_of_haplotypes_plus_a_bit} /lizardfs/guarracino/chromosome_communities/graphs";


# Full acrocentric pangenome graph
sbatch -p workers -c 48 --job-name acropggb --wrap "hostname; cd /scratch && pggb -i /lizardfs/guarracino/chromosome_communities/pq_contigs/s5k/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz -o chrACRO+refs.pq_contigs_s5k.1kbps.hg002prox.hg002hifi.s25k.p90.n300 -t 48 -s 25k -p 90 -F 0.001 -n 300 -k 311 -G 13117,13219 -O 0.03 -T 48 -V chm13:#,grch38:#; mv /scratch/chrACRO+refs.pq_contigs_s5k.1kbps.hg002prox.hg002hifi.s25k.p90.n300 /lizardfs/guarracino/chromosome_communities/graphs";
```


### Untangling

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/untangle

path_targets_txt=/lizardfs/guarracino/chromosome_communities/untangle/chm13.target_paths.txt
grep chm13 /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.fa.gz.fai | cut -f 1 > $path_targets_txt

# Graph
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.og
prefix=$(basename $path_input_og .og)

RUN_ODGI=/home/guarracino/tools/odgi/bin/odgi-454197fa29b772050c3135d5de47c816ce38e62c
```


Untangle with respect to all acrocentric chromosomes and emit the cut points:

```shell
for e in 50000; do
  for m in 1000; do
    path_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.j0.n100.bed.gz
    path_cut_points_txt=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.j0.n100.cut_points.txt
    
    if [[ ! -s ${path_cut_points_txt} ]]; then
      echo "-e $e -m $m"
      sbatch -p workers -c 48 --job-name acrountangle --wrap "\time -v $RUN_ODGI untangle -t 48 -P -i $path_input_og -R $path_targets_txt -e $e -m $m --cut-points-output $path_cut_points_txt -j 0 -n 100 | pigz -c > $path_bed_gz"
    fi;
  done
done
```


Fix best hits (if there are multiple best hits, put as first the target-chromosome of origin of the contig):

```shell
for e in 50000; do
  for m in 1000 ; do
    path_ref_fixed_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.j0.n100.fixed.bed.gz
    if [[ ! -s ${path_ref_fixed_bed_gz} ]]; then
      echo "-e $e -m $m"
      
      path_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.j0.n100.bed.gz
      python3 /lizardfs/guarracino/chromosome_communities/scripts/fix_best_hit.py \
        $path_bed_gz \
        <(cat \
            <(cat assemblies/partitioning/*.partitions.tsv | sed 's/chr/chm13#chr/') \
            <(cut -f 1,6  /lizardfs/guarracino/chromosome_communities/assemblies/partitioning_pq/HG002.*.vs.ref.p90.N.paf) \
            <(echo -e chm13#chr13"\t"chm13#chr13) \
            <(echo -e chm13#chr14"\t"chm13#chr14) \
            <(echo -e chm13#chr15"\t"chm13#chr15) \
            <(echo -e chm13#chr21"\t"chm13#chr21) \
            <(echo -e chm13#chr22"\t"chm13#chr22) \
            <(echo -e grch38#chr13"\t"chm13#chr13) \
            <(echo -e grch38#chr14"\t"chm13#chr14) \
            <(echo -e grch38#chr15"\t"chm13#chr15) \
            <(echo -e grch38#chr21"\t"chm13#chr21) \
            <(echo -e grch38#chr22"\t"chm13#chr22) \
            <(echo -e HG002#MAT#chr13.prox"\t"chm13#chr13) \
            <(echo -e HG002#MAT#chr14.prox"\t"chm13#chr14) \
            <(echo -e HG002#MAT#chr15.prox"\t"chm13#chr15) \
            <(echo -e HG002#MAT#chr21.prox"\t"chm13#chr21) \
            <(echo -e HG002#MAT#chr22.prox"\t"chm13#chr22) \
            <(echo -e HG002#PAT#chr13.prox"\t"chm13#chr13) \
            <(echo -e HG002#PAT#chr14.prox"\t"chm13#chr14) \
            <(echo -e HG002#PAT#chr15.prox"\t"chm13#chr15) \
            <(echo -e HG002#PAT#chr21.prox"\t"chm13#chr21) \
            <(echo -e HG002#PAT#chr22.prox"\t"chm13#chr22)) | tr ' ' '\t' | pigz -c -9 > $path_ref_fixed_bed_gz
    fi;
  done
done

# Fixed 7661 hits covering 11171215 bps on the queries.
```


Untangle with respect to a single acrocentric chromosome, using always the same cut points:

```shell
for e in 50000; do
  for m in 1000; do
    echo "-e $e -m $m"
      
    cat $path_targets_txt | while read ref; do
      echo $ref
      
      path_ref_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$ref.e$e.m$m.j0.n100.bed.gz
      if [[ ! -s ${path_ref_bed_gz} ]]; then
        path_cut_points_txt=/lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.j0.n100.cut_points.txt
        
        sbatch -p workers -c 48 --job-name acrountangle --wrap "\time -v $RUN_ODGI untangle -t 48 -P -i $path_input_og -r $ref -e $e -m $m --cut-points-input $path_cut_points_txt -j 0 -n 100 | pigz -c > $path_ref_bed_gz"
      fi;
    done
  done
done
```


Grounding (applying filters) and annotation:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/untangle/grounded

for e in 50000; do
  for m in 1000 ; do
    cat $path_targets_txt | while read ref; do            
      path_grounded_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.grounded.tsv.gz
            
      if [[ ! -s ${path_grounded_tsv_gz} ]]; then
        echo "-e $e -m $m $ref grounding"

        # Grounding
        ( echo query query.begin query.end target target.begin target.end jaccard strand self.coverage nth.best ref ref.begin ref.end ref.jaccard ref.nth.best | tr ' ' '\t'
          join \
            <(zcat /lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.chm13#ACRO.e$e.m$m.j0.n100.fixed.bed.gz | awk -v j=0 -v n=5 '{ if($7 >= j && $10 <= n) {print $1"_"$2, $0} }' | tr ' ' '\t' | sort -k 1,1) \
            <(zcat /lizardfs/guarracino/chromosome_communities/untangle/$prefix.untangle.$ref.e$e.m$m.j0.n100.bed.gz | awk -v j=0 -v n=10 '{ if($7 >= j && $10 <= n) {print $1"_"$2, $0} }' | tr ' ' '\t' | sort -k 1,1) | \
          tr ' ' '\t' | grep -v '^#' | cut -f 2- | cut -f -10,14-17,20 | sort -k 1,3 -k 7,7nr -k 10,10n -k 14,14nr -k 15,15n ) | tr ' ' '\t' | pigz -c -9 > ${path_grounded_tsv_gz}
      fi;

      # Take pq-untangling contigs
               path_grounded_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.grounded.tsv.gz
      path_grounded_pq_touching_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.grounded.pq_touching.tsv
          
      if [[ ! -s ${path_grounded_pq_touching_tsv}.gz ]]; then
        echo "-e $e -m $m $ref filtering&annotation"

        # "p-touching contigs" intersected "q-touching contigs"
        comm -12 \
          <(bedtools intersect \
            -a <(zcat ${path_grounded_tsv_gz} | awk -v OFS="\t" '{print $11,$12,$13,$1, "", "+"}' | sed '1d' | bedtools sort) \
            -b <(grep $ref /lizardfs/guarracino/chromosome_communities/assemblies/partitioning_pq/p_arms.bed | bedtools sort) | \
            #awk '$3-$2+1>=1000' | \
            cut -f 4 | sort | uniq) \
          <(bedtools intersect \
            -a <(zcat ${path_grounded_tsv_gz} | awk -v OFS="\t" '{print $11,$12,$13,$1, "", "+"}' | sed '1d' | bedtools sort) \
            -b <(grep $ref /lizardfs/guarracino/chromosome_communities/assemblies/partitioning_pq/q_arms.bed | bedtools sort) | \
            #awk '$3-$2+1>=1000' | \
            cut -f 4 | sort | uniq) | \
            grep 'chm13#\|grch38#' -v > $ref.tmp.txt
        #zcat ${path_grounded_tsv_gz} | cut -f 1 | sed '1d' | sort | uniq >  $ref.tmp.txt # ALL
      
        # Consider only chromosome-partitioned contigs
        comm -12 \
          <(sort $ref.tmp.txt) \
          <(cut -f 1 /lizardfs/guarracino/chromosome_communities/pq_contigs/s5k/$ref.pq_contigs.1kbps.hg002prox.fa.gz.fai | sort)  \
          > $ref.tmp2.txt
        rm $ref.tmp.txt
        
        ##########################################################################################
        # To keep the short (bad) HiFi-based HG002 contigs
        chr=$( echo $ref | sed 's/chm13#//g')
        # To avoid possible duplicates (for HG002's pq-contigs)
        comm -13 \
          <(cat $ref.tmp2.txt | sort) \
          <(cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning_pq/HG002.*.vs.ref.p90.N.paf | grep $chr | awk '$2 >= 300000' | cut -f 1 | sort) | \
          sort | uniq > $ref.tmp3.txt
        cat $ref.tmp3.txt >> $ref.tmp2.txt
        rm $ref.tmp3.txt
        
#            # In case we have to fix a partitioning
#            if [ $ref == "chm13#chr13" ]; then
##              echo "Remove wrongly partitioned contig: HG002#1#h1tg000013l"
##              grep 'HG002#1#h1tg000013l' -v $ref.tmp2.txt > $ref.tmp3.txt
##              echo "Add wrongly partitioned contig: HG002#1#h1tg000260l"
##              echo 'HG002#1#h1tg000260l' >> $ref.tmp3.txt
#
#              echo "Remove wrongly partitioned contig: HG002#1#JAHKSE010000013.1"
#              grep 'HG002#1#JAHKSE010000013.1' -v $ref.tmp2.txt > $ref.tmp3.txt
#              echo "Add wrongly partitioned contig: HG002#1#JAHKSE010000214.1"
#              echo 'HG002#1#JAHKSE010000214.1' >> $ref.tmp3.txt             
#              
#              rm $ref.tmp2.txt && mv $ref.tmp3.txt $ref.tmp2.txt
#            fi           
        ##########################################################################################
      
        #### Put header (with a new 'target' column), take intersection, and re-add other acrocentric references
        cat \
          <(zcat ${path_grounded_tsv_gz} | head -n 1 | awk -v OFS='\t' '{print $0, "grounded.target"}') \
          <(zgrep ${path_grounded_tsv_gz} -f $ref.tmp2.txt | awk -v OFS='\t' -v ref=$ref '{print $0,ref}') \
          <(zgrep '^grch38\|^chm13' ${path_grounded_tsv_gz}  | awk -v OFS='\t' -v ref=$ref '{print $0,ref}') \
            > ${path_grounded_pq_touching_tsv}
        rm $ref.tmp2.txt
              
        # Add annotation
        zgrep rDNA /lizardfs/guarracino/chromosome_communities/data/chm13.CentromeresAndTelomeres.CenSatAnnotation.txt.gz | \
          sed 's/chr/chm13#chr/g' | \
          grep $ref | \
          awk -v OFS='\t' -v ref=$ref '{print $1"#"$4,".",".",$1,".",".","1","+","1","1",ref,$2,$3,"1","1",ref}' >> ${path_grounded_pq_touching_tsv}
        cat /lizardfs/guarracino/chromosome_communities/data/chm13.centromeres.approximate.bed | \
          bedtools merge | \
          grep 'chr13\|chr14\|chr15\|chr21\|chr22' | \
          sed 's/chr/chm13#chr/g'  | \
          grep $ref | \
          awk -v OFS='\t' -v ref=$ref '{print $1"#centromere",".",".",$1,".",".","1","+","1","1",ref,$2,$3,"1","1",ref}' >> ${path_grounded_pq_touching_tsv}
              
        pigz ${path_grounded_pq_touching_tsv} -f -9
      fi;
    done
  done
done
```


Remove unreliable regions:

```shell
rm x.tsv # Cleaning

for e in 50000; do
  for m in 1000; do
    cat $path_targets_txt | while read ref; do
               path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.grounded.pq_touching.tsv.gz
      path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.grounded.pq_touching.reliable.tsv.gz

      if [[ ! -s ${path_grounded_pq_touching_reliable_tsv_gz} ]]; then
        echo "-e $e -m $m $ref"
        
        # Skip verkko's and bakeoff's contigs because we don't have the unreliable regions for them
        zcat $path_grounded_pq_touching_tsv_gz | sed '1d' | grep 'HG002#MAT\|HG002#PAT\|HG01978#MAT\|HG01978#PAT\|bakeoff' -v | cut -f 1 | sort | uniq | while read CONTIG; do
          SAMPLE=$( echo $CONTIG | cut -f 1 -d '#')

          path_unreliable_bed=/lizardfs/guarracino/HPRC/annotations/unreliable/$SAMPLE.hifi.flagger_final.simplified.unreliable_only.bed
          if [[ -s $path_unreliable_bed ]]; then
            #echo $CONTIG "--->" $SAMPLE
            
            zgrep "^$CONTIG" $path_grounded_pq_touching_tsv_gz | sort -k 2n | awk -v OFS='\t' '{print($1,$2,$3,$4"_"$5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11"_"$12"_"$13"_"$14"_"$15"_"$16)}' > x.bed
            grep $CONTIG $path_unreliable_bed > y.bed
            # -A: remove entire feature if any overlap
            bedtools subtract -a x.bed -b y.bed -A |\
              awk -v OFS='\t' '{split($4, a, "_"); print($1,$2,$3,a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[11],a[12],a[13])}' >> x.tsv
            rm x.bed y.bed
          else
            zgrep "^$CONTIG" $path_grounded_pq_touching_tsv_gz | sort -k 2n >> x.tsv
          fi
        done
        
        # Re-take verkko's and bakeoff's untangled regions
        zcat $path_grounded_pq_touching_tsv_gz | sed '1d' | grep 'HG002#MAT\|HG002#PAT\|HG01978#MAT\|HG01978#PAT\|bakeoff' >> x.tsv
        
        cat \
          <( zcat $path_grounded_pq_touching_tsv_gz | head -n 1 ) \
          <( sort -k 1,3 -k 7,7nr -k 10,10n -k 14,14nr -k 15,15n x.tsv) | pigz -c -9 > $path_grounded_pq_touching_reliable_tsv_gz
        rm x.tsv
      fi;
    done
  done
done
```


Annotated plots:

```shell
# Dependencies
#guix install r
#guix install r-ggplot2
#guix install r-ggforce
#guix install r-tidyverse
#guix install r-png
#guix install r-ggpubr
#guix install r-scales
#guix install poppler # For pdfunite

# Annotated plots with the first hit
for e in 50000; do
  for m in 1000; do
    for refn in 1 10; do
      (seq 13 15; seq 21 22) | while read i; do
        echo "-e $e -m $m -refn $refn chr$i"

        path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chr${i}.e$e.m$m.grounded.pq_touching.reliable.tsv.gz
        PREFIX=$(basename $path_grounded_pq_touching_reliable_tsv_gz .tsv.gz);
        
        Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation.R \
          $path_grounded_pq_touching_reliable_tsv_gz \
          0 25000000 \
          90 0.7 \
          0 \
          1 $refn \
          $i \
          0.9 \
          <(zgrep '^HG002#1\|^HG002#2' -v $path_grounded_pq_touching_reliable_tsv_gz | sed '1d' | cut -f 1 | sort | uniq) \
          /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \
          /lizardfs/guarracino/chromosome_communities/untangle/grounded/$PREFIX.n1.nref${refn}.pdf
      done
      
      # Merge chromosomes's PDF files
      /gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite \
        /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chr*.e$e.m$m.grounded.pq_touching.reliable.n1.nref${refn}.pdf \
        /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.n1.nref${refn}.merged.pdf
      rm /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chr*.e$e.m$m.grounded.pq_touching.reliable.n1.nref${refn}.pdf
    done
  done
done

# Annotated plots with the first 5 hits
for e in 50000; do
  for m in 1000; do
    for refn in 1 10; do
      (seq 13 15; seq 21 22) | while read i; do
        echo "-e $e -m $m -refn $refn chr$i"
        
        path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chr${i}.e$e.m$m.grounded.pq_touching.reliable.tsv.gz
        PREFIX=$(basename $path_grounded_pq_touching_reliable_tsv_gz .tsv.gz);
        
        Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation.R \
          $path_grounded_pq_touching_reliable_tsv_gz \
          0 25000000 \
          91 0.8 \
          0.8 \
          5 $refn \
          $i \
          0.9 \
          <(zgrep '^HG002#1\|^HG002#2' -v $path_grounded_pq_touching_reliable_tsv_gz | sed '1d' | cut -f 1 | sort | uniq) \
          /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \
          /lizardfs/guarracino/chromosome_communities/untangle/grounded/$PREFIX.n5.nref${refn}.pdf
      done
      
      # Merge chromosomes's PDF files
      /gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite \
        /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chr*.e$e.m$m.grounded.pq_touching.reliable.n5.nref${refn}.pdf \
        /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.n5.nref${refn}.merged.pdf
      rm /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chr*.e$e.m$m.grounded.pq_touching.reliable.n5.nref${refn}.pdf
    done
  done
done
```


Compute support by considering HiFi-only contigs anchored to the q-arms (so no HG002-HiFi-only) and HG002-verkko:

```shell
# Merge files for all acrocentric chromosomes (used for computing the support and the histogram length)
for e in 50000; do
  for m in 1000; do
    path_grounded_pq_touching_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.ALL.e$e.m$m.grounded.pq_touching.reliable.tsv.gz
    if [[ ! -s ${path_grounded_pq_touching_reliable_all_chromosomes_tsv_gz} ]]; then
      cat \
        <(zcat /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chr*.e$e.m$m.grounded.pq_touching.reliable.tsv.gz | grep 'self.coverage' -m 1) \
        <(zcat /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chr*.e$e.m$m.grounded.pq_touching.reliable.tsv.gz | grep 'self.coverage' -v ) |\
        pigz -c -9 > $path_grounded_pq_touching_reliable_ALL_tsv_gz
    fi;
  done
done

#Take acrocentric chromosome lengths
grep '^chm13' /lizardfs/guarracino/chromosome_communities/assemblies/chrA.pan+HG002chrAprox.fa.gz.fai | cut -f 1,2 \
  > /lizardfs/guarracino/chromosome_communities/chm13#ACRO.len.tsv

# Support
# guix install r-ggridges
for e in 50000; do
  for m in 1000; do
    for eid in 0.900 0.950 0.975 0.995 1.000; do
      eid_str=$(echo $eid | sed 's/\.//g')
      for refn in 1; do
        echo "-e $e -m $m $eid -refn $refn"

        path_grounded_pq_touching_reliable_ALL_support_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.support.eid${eid_str}.n1.nref${refn}.tsv.gz
        if [[ ! -s  $path_grounded_pq_touching_reliable_ALL_support_tsv_gz ]]; then
          path_grounded_pq_touching_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.ALL.e$e.m$m.grounded.pq_touching.reliable.tsv.gz
          python3 /lizardfs/guarracino/chromosome_communities/scripts/support.py \
            $path_grounded_pq_touching_reliable_ALL_tsv_gz \
            /lizardfs/guarracino/chromosome_communities/chm13#ACRO.len.tsv 1 $refn \
            $eid \
            <( zgrep '^chm13\|^grch38\|^HG002#1\|HG002#2\|^HG01978#MAT\|^HG01978#PAT\|bakeoff' $path_grounded_pq_touching_reliable_ALL_tsv_gz -v | sed '1d' | cut -f 1 | sort | uniq ) | \
            pigz -c -9 > $path_grounded_pq_touching_reliable_ALL_support_tsv_gz
        fi
    
        path_grounded_pq_touching_reliable_ALL_support_dedup_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.support.dedup.eid${eid_str}.n1.nref${refn}.tsv.gz
        if [[ ! -s $path_grounded_pq_touching_reliable_ALL_support_dedup_tsv_gz ]]; then
          python3 /lizardfs/guarracino/chromosome_communities/scripts/support_deduplication.py \
            $path_grounded_pq_touching_reliable_ALL_support_tsv_gz | \
            pigz -c > $path_grounded_pq_touching_reliable_ALL_support_dedup_tsv_gz
        fi
        
        PREFIX=$(basename $path_grounded_pq_touching_reliable_ALL_support_dedup_tsv_gz .tsv.gz);
        (seq 13 15; seq 21 22) | while read i; do
          echo chr$i
          Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_collapsed_with_annotation.R \
            $path_grounded_pq_touching_reliable_ALL_support_dedup_tsv_gz \
            0 25000000 \
            84 \
            $i \
            /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \
            /lizardfs/guarracino/chromosome_communities/untangle/grounded/$PREFIX.chr$i
        done
        
        # Merge chromosomes's PDF files
        /gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite \
          /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.support.dedup.eid${eid_str}.n1.nref${refn}.chr*.separated.pdf \
          /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.support.dedup.eid${eid_str}.n1.nref${refn}.separated.merged.pdf
        rm /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.support.dedup.eid${eid_str}.n1.nref${refn}.chr*.separated.pdf
        /gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite \
          /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.support.dedup.eid${eid_str}.n1.nref${refn}.chr*.together.pdf \
          /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.support.dedup.eid${eid_str}.n1.nref${refn}.together.merged.pdf
        rm /lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.support.dedup.eid${eid_str}.n1.nref${refn}.chr*.together.pdf 
      done
    done
  done
done
```


Statistics on removed regions (available only for HiFi-only samples):

```shell
path_grounded_pq_touching_reliable_stats_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.grounded.pq_touching.reliable.stats.tsv
    
echo -e "ground\te\tm\testimated.identity\tn\tref.n\tcontig\tuntangled.size\treliable.untangled.size\tfraction.removed" > $path_grounded_pq_touching_reliable_stats_tsv

n=1
refn=1    
for e in 50000; do
  for m in 1000; do
    for eid in 0.900 0.950 0.975 0.995 1.000; do     
      cat $path_targets_txt | while read ref; do
        echo "-e $e -m $m $eid $ref"
                 path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.grounded.pq_touching.tsv.gz
        path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.grounded.pq_touching.reliable.tsv.gz

        zgrep '^chm13\|^grch38\|^HG002#MAT\|^HG002#PAT\|^HG01978#MAT\|^HG01978#PAT\|bakeoff' -v $path_grounded_pq_touching_tsv_gz | sed '1d' | cut -f 1 | sort | uniq | while read CONTIG; do
          UNTANGLED_SIZE=$( zgrep "^$CONTIG" $path_grounded_pq_touching_tsv_gz | awk -v eid=$eid -v n=$n -v refn=$refn 'current_eid=exp((1.0 + log(2.0 * $7/(1.0+$7)))-1.0); current_eid >= eid && $10 == n && $15 == refn' | cut -f 1,2,3 | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' )
          UNTANGLED_SIZE_RELIABLE=$( zgrep "^$CONTIG" $path_grounded_pq_touching_reliable_tsv_gz | awk -v eid=$eid -v n=$n -v refn=$refn 'current_eid=exp((1.0 + log(2.0 * $7/(1.0+$7)))-1.0); current_eid >= eid && $10 == n && $15 == refn' | cut -f 1,2,3 | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' )
        
          FRACTION_REMOVED=$(echo "scale=4; 1 - $UNTANGLED_SIZE_RELIABLE/$UNTANGLED_SIZE" | bc)
      
          echo $ref $e $m $eid $n $refn $CONTIG $UNTANGLED_SIZE $UNTANGLED_SIZE_RELIABLE $FRACTION_REMOVED | tr ' ' '\t' >> $path_grounded_pq_touching_reliable_stats_tsv
        done
      done
    done
  done
done

# Check total over the contigs
(echo "untangle.space.bp untangle.space.reliable.bp untangle.space.unreliable.bp, fraction.untangle.space.reliable" ; \
  awk '$2 == 50000 && $3 == 1000 && $4 == 0.900 && $5 == 1 && $6 == 1' $path_grounded_pq_touching_reliable_stats_tsv | \
    awk -v OFS='\t' -F'\t' 'BEGIN{UNTANGLED_SIZE=0; UNTANGLED_SIZE_RELIABLE=0}{ UNTANGLED_SIZE+=$8; UNTANGLED_SIZE_RELIABLE+=$9 }END{print UNTANGLED_SIZE,UNTANGLED_SIZE_RELIABLE,UNTANGLED_SIZE-UNTANGLED_SIZE_RELIABLE, UNTANGLED_SIZE_RELIABLE/UNTANGLED_SIZE}' )

(head -n 1  $path_grounded_pq_touching_reliable_stats_tsv; \
  awk '$2 == 50000 && $3 == 1000 && $4 == 0.900 && $5 == 1 && $6 == 1' $path_grounded_pq_touching_reliable_stats_tsv) | \
  cut -f 7,8,9,10 > SuppTable2.eid0900.tsv
```


Statistics on untangled segment lengths by considering HiFi-only contigs anchored to the q-arms (so no HG002-HiFi-only) and HG002-verkko:

```shell
for e in 50000; do
  for m in 1000; do
    path_grounded_pq_touching_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.ALL.e$e.m$m.grounded.pq_touching.reliable.tsv.gz
    PREFIX=$(basename $path_grounded_pq_touching_reliable_ALL_tsv_gz .tsv.gz);
      
    for eid in 0.900 0.950 0.975 0.995 1.000; do
      eid_str=$(echo $eid | sed 's/\.//g')
    
      Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangled_segment_histogram.R \
        $path_grounded_pq_touching_reliable_ALL_tsv_gz \
        0 25000000 \
        60 15 \
        $(echo "$e + 15000" | bc) \
        1 1 \
        $eid \
        <( zgrep '^chm13\|^grch38\|^HG002#1\|HG002#2\|^HG01978#MAT\|^HG01978#PAT\|bakeoff' $path_grounded_pq_touching_reliable_ALL_tsv_gz -v | sed '1d' | cut -f 1 | sort | uniq ) \
        /lizardfs/guarracino/chromosome_communities/untangle/grounded/$PREFIX.eid${eid_str}.n1.nref1.histogram.pdf
    done
  done
done
```


Estimate regions that can recombine using multi-hit untangled regions:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/untangle/grounded/recombinant_regions/

# Identify regions with multiple good enough hits
for e in 50000; do
  for m in 1000; do
    path_grounded_pq_touching_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.ALL.e$e.m$m.grounded.pq_touching.reliable.tsv.gz
    PREFIX=$(basename $path_grounded_pq_touching_reliable_ALL_tsv_gz .tsv.gz)

    zgrep '^chm13\|^grch38\|^HG002#1\|HG002#2\|^HG01978#MAT\|^HG01978#PAT\|bakeoff' $path_grounded_pq_touching_reliable_ALL_tsv_gz -v | sed '1d' | cut -f 1 | sort | uniq \
      > /lizardfs/guarracino/chromosome_communities/untangle/grounded/recombinant_regions/$PREFIX.query_to_consider.txt
   
    for sc in 0 1.1; do
      for eid in 0.900 0.950 0.975 0.995 1.000; do
        eid_str=$(echo $eid | sed 's/\.//g')
        sc_str=$(echo $sc | sed 's/\.//g')
        echo $e $m $sc $eid
          
        path_recombinant_regions_bed=/lizardfs/guarracino/chromosome_communities/untangle/grounded/recombinant_regions/$PREFIX.recombinant_regions.sc${sc_str}.eid${eid_str}.bed
        if [[ ! -s $path_recombinant_regions_bed ]]; then
          python3 /lizardfs/guarracino/chromosome_communities/scripts/recombination_proxy_ranges.py \
            $path_grounded_pq_touching_reliable_ALL_tsv_gz \
            $eid $sc \
            /lizardfs/guarracino/chromosome_communities/untangle/grounded/recombinant_regions/$PREFIX.query_to_consider.txt \
            > $path_recombinant_regions_bed
        fi
      done
    done
  done
done

# Collect values in grounded reference space (TO PUT IN A FILE AND RUN A JOB FOR IT)
########################################################################################################################
#!/bin/sh

cd /scratch
awk -v OFS='\t' '{print($1,"0",$2)}' /lizardfs/guarracino/chromosome_communities/chm13#ACRO.len.tsv > chm13.bed

path_targets_txt=/lizardfs/guarracino/chromosome_communities/untangle/chm13.target_paths.txt
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.og
prefix=$(basename $path_input_og .og)

for e in 50000; do
  for m in 1000; do
    path_grounded_pq_touching_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.ALL.e$e.m$m.grounded.pq_touching.reliable.tsv.gz
    PREFIX=$(basename $path_grounded_pq_touching_reliable_ALL_tsv_gz .tsv.gz)

    path_recombinant_regions_table_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/recombinant_regions/$PREFIX.recombinant_regions.table.tsv
    path_recombinant_regions_table_sizes_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/recombinant_regions/$PREFIX.recombinant_regions.table.sizes.tsv
    path_recombinant_regions_table_with_counts_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/recombinant_regions/$PREFIX.recombinant_regions.table.counts.tsv
    
    rm rg.txt
    rm $path_recombinant_regions_table_tsv $path_recombinant_regions_table_sizes_tsv $path_recombinant_regions_table_with_counts_tsv
    for sc in 0 1.1; do
      for eid in 0.900 0.950 0.975 0.995 1.000; do
        eid_str=$(echo $eid | sed 's/\.//g')
        sc_str=$(echo $sc | sed 's/\.//g')
        echo $e $m $sc $eid
          
        path_recombinant_regions_bed=/lizardfs/guarracino/chromosome_communities/untangle/grounded/recombinant_regions/$PREFIX.recombinant_regions.sc${sc_str}.eid${eid_str}.bed

        bedtools merge -i <(cut -f 4,5,6 $path_recombinant_regions_bed | sed '1d' | bedtools sort ) | \
          awk -v sc=$sc -v eid=$eid -v OFS='\t' '{print(sc,eid,$1,$2,$3)}' >> $path_recombinant_regions_table_tsv
          
        bedtools merge -i <(cut -f 4,5,6 $path_recombinant_regions_bed | sed '1d' | bedtools sort ) | \
          awk -v sc=$sc -v eid=$eid -v OFS='\t' '{SUM+=$3-$2}END{print(sc,eid,SUM)}' >> $path_recombinant_regions_table_sizes_tsv

        # For each sample, merge intervals with respect to the grounded reference
        sed '1d' $path_recombinant_regions_bed | cut -f 1 | sort | uniq | while read CONTIG; do
          grep "^$CONTIG" $path_recombinant_regions_bed | cut -f 4,5,6 | bedtools sort | bedtools merge >> rg.txt
        done
        
        # For each grounded reference position, count how many sample support it
        # -d: Report the depth at each position in each A feature.
        # The awk script is to get intervals where grounded reference and counts is constant.
        bedtools coverage -a chm13.bed -b rg.txt -d | \
          python3 /lizardfs/guarracino/chromosome_communities/scripts/compress_coverage_info.py | \
          awk -v sc=$sc -v eid=$eid -v OFS='\t' '{print(sc,eid,$1,$2,$3,$4)}' \
          >> $path_recombinant_regions_table_with_counts_tsv
        rm rg.txt
      done
    done
  done
done

rm chm13.bed
########################################################################################################################

# Plots
for e in 50000; do
  for m in 1000; do
    path_grounded_pq_touching_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.ALL.e$e.m$m.grounded.pq_touching.reliable.tsv.gz
    PREFIX=$(basename $path_grounded_pq_touching_reliable_ALL_tsv_gz .tsv.gz)
    
    path_recombinant_regions_table_with_counts_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/recombinant_regions/$PREFIX.recombinant_regions.table.counts.tsv
    
    (seq 13 15; seq 21 22) | while read i; do
      echo "-e $e -m $m chr$i"
        
      Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_recombinant_regions_with_annotation.R \
        $path_recombinant_regions_table_with_counts_tsv \
        0 25000000 \
        90 15 \
        $i \
        /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \
        /lizardfs/guarracino/chromosome_communities/untangle/grounded/recombinant_regions/$PREFIX.recombinant_regions.table.chr${i}.pdf
    done
    
    # Merge chromosomes's PDF files
    /gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite \
      /lizardfs/guarracino/chromosome_communities/untangle/grounded/recombinant_regions/$PREFIX.recombinant_regions.table.chr*.pdf \
      /lizardfs/guarracino/chromosome_communities/untangle/grounded/recombinant_regions/$PREFIX.recombinant_regions.table.merged.pdf
    rm /lizardfs/guarracino/chromosome_communities/untangle/grounded/recombinant_regions/$PREFIX.recombinant_regions.table.chr*.pdf
  done
done

for e in 50000; do
  for m in 1000; do
    path_grounded_pq_touching_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.ALL.e$e.m$m.grounded.pq_touching.reliable.tsv.gz
    PREFIX=$(basename $path_grounded_pq_touching_reliable_ALL_tsv_gz .tsv.gz)
    
    path_recombinant_regions_table_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/recombinant_regions/$PREFIX.recombinant_regions.table.tsv
    path_recombinant_regions_table_sizes_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/recombinant_regions/$PREFIX.recombinant_regions.table.sizes.tsv
    path_recombinant_regions_table_with_counts_tsv=/lizardfs/guarracino/chromosome_communities/untangle/grounded/recombinant_regions/$PREFIX.recombinant_regions.table.counts.tsv
    
    Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_threshold_vs_recombinant.R \
      $path_recombinant_regions_table_sizes_tsv \
      /lizardfs/guarracino/chromosome_communities/untangle/grounded/recombinant_regions/threshold_vs_recombinant_regions.pdf
  done
done




# Unmerged query segments
cat x.bed | wc -l

# Merged query segments: take query columns, merge contiguos intervals and count them
bedtools merge -i <(cut -f 1,2,3 x.bed | sed '1d' ) | wc -l

# Merged ground segments
bedtools merge -i <(cut -f 4,5,6 x.bed | sed '1d' | bedtools sort ) | wc -l
```


How many pq-contigs appear to cross the rDNA array (touches both sides)?

```shell
n=1
refn=1    
for e in 50000; do
  for m in 1000; do
    path_grounded_pq_touching_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.ALL.e$e.m$m.grounded.pq_touching.reliable.tsv.gz
    PREFIX=$(basename $path_grounded_pq_touching_reliable_ALL_tsv_gz .tsv.gz)
    
    for eid in 0.900 0.950 0.975 0.995 1.000; do
      eid_str=$(echo $eid | sed 's/\.//g')
      
      echo "-e $e -m $m $eid"
      
      comm -12 \
        <( bedtools intersect \
            -a <( zgrep '^chm13\|^grch38\|^HG002#1\|HG002#2\|^HG01978#MAT\|^HG01978#PAT\|bakeoff' -v $path_grounded_pq_touching_reliable_ALL_tsv_gz | sed '1d' | \
                    awk -v eid=$eid -v n=$n -v refn=$refn 'current_eid=exp((1.0 + log(2.0 * $7/(1.0+$7)))-1.0); current_eid >= eid && $10 == n && $15 == refn' | \
                    awk -v OFS='\t' '{print($11,$12,$13,$1)}' | bedtools sort ) \
            -b <( zgrep rDNA /lizardfs/guarracino/chromosome_communities/data/chm13.CentromeresAndTelomeres.CenSatAnnotation.txt.gz | \
                  awk -v OFS='\t' '{print("chm13#"$1,"0",$2)}' | bedtools sort ) | cut -f 4 | sort | uniq)  \
        <( bedtools intersect \
            -a <( zgrep '^chm13\|^grch38\|^HG002#1\|HG002#2\|^HG01978#MAT\|^HG01978#PAT\|bakeoff' -v $path_grounded_pq_touching_reliable_ALL_tsv_gz | sed '1d' | \
                    awk -v eid=$eid -v n=$n -v refn=$refn 'current_eid=exp((1.0 + log(2.0 * $7/(1.0+$7)))-1.0); current_eid >= eid && $10 == n && $15 == refn' | \
                    awk -v OFS='\t' '{print($11,$12,$13,$1)}' | bedtools sort ) \
            -b <( bedtools complement \
                    -i <( zgrep rDNA /lizardfs/guarracino/chromosome_communities/data/chm13.CentromeresAndTelomeres.CenSatAnnotation.txt.gz | \
                              awk -v OFS='\t' '{print("chm13#"$1,"0",$2)}' | bedtools sort ) \
                    -g <( grep 'chr13\|chr14\|chr15\|chr21\|chr22' /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa.fai | cut -f 1,2 | sort ) ) | cut -f 4 | sort | uniq ) | \
        sort | uniq > /lizardfs/guarracino/chromosome_communities/untangle/grounded/$PREFIX.eid${eid_str}.rDNA_crossing_contigs.tsv
    done
  done
done
```

Self coverage plots:

```shell
for e in 50000; do
  for m in 1000; do
    path_grounded_pq_touching_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.ALL.e$e.m$m.grounded.pq_touching.reliable.tsv.gz
    PREFIX=$(basename $path_grounded_pq_touching_reliable_ALL_tsv_gz .tsv.gz);
      
    for eid in 0.900; do
      eid_str=$(echo $eid | sed 's/\.//g')
    
      (seq 13 15; seq 21 22) | while read i; do
        Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_self_coverage_with_annotation.R \
          $path_grounded_pq_touching_reliable_ALL_tsv_gz \
          0 25000000 \
          90 4 0.6 \
          1 1 \
          $i \
          $eid \
          <( zgrep '^chm13\|^grch38\|^HG002#1\|HG002#2\|^HG01978#MAT\|^HG01978#PAT\|bakeoff' $path_grounded_pq_touching_reliable_ALL_tsv_gz -v | sed '1d' | cut -f 1 | sort | uniq ) \
          /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \
          /lizardfs/guarracino/chromosome_communities/untangle/grounded/$PREFIX.eid${eid_str}.n1.nref1.self_coverage.chr${i}.pdf
      done
    done
  done
done
```

[//]: # (Plot with manually selected paths:)
[//]: # ()
[//]: # (```shell)
[//]: # (e=50000)
[//]: # (m=1000)
[//]: # ()
[//]: # (i=13)
[//]: # (path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chr${i}.e$e.m$m.grounded.pq_touching.reliable.tsv.gz)
[//]: # (PREFIX=$&#40;basename $path_grounded_pq_touching_reliable_tsv_gz .tsv.gz&#41;;)
[//]: # (Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation.R \)
[//]: # (  $path_grounded_pq_touching_reliable_tsv_gz \)
[//]: # (  0 25000000 \)
[//]: # (  91 0.8 \)
[//]: # (  0 \)
[//]: # (  1 1 \)
[//]: # (  $i \)
[//]: # (  0.9 \)
[//]: # (  <&#40;echo chm13#chr13 grch38#chr13 HG002#MAT#chr13.prox HG002#PAT#chr13.prox HG01361#2#JAGYYW010000010.1 HG01978#1#JAGYVS010000056.1 HG02486#1#JAGYVM010000043.1 HG03540#2#JAGYVX010000153.1 | tr ' ' '\n'&#41; \)
[//]: # (  /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \)
[//]: # (  ~/$PREFIX.n1.subset.pdf)
[//]: # (Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation.R \)
[//]: # (  $path_grounded_pq_touching_reliable_tsv_gz \)
[//]: # (  0 25000000 \)
[//]: # (  91 0.8 \)
[//]: # (  0.8 \)
[//]: # (  5 1 \)
[//]: # (  $i \)
[//]: # (  0.9 \)
[//]: # (  <&#40;echo chm13#chr13 grch38#chr13 HG002#MAT#chr13.prox HG002#PAT#chr13.prox HG01361#2#JAGYYW010000010.1 HG01978#1#JAGYVS010000056.1 HG02486#1#JAGYVM010000043.1 HG03540#2#JAGYVX010000153.1 | tr ' ' '\n'&#41; \)
[//]: # (  /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \)
[//]: # (  /lizardfs/guarracino/chromosome_communities/untangle/grounded/$PREFIX.n5.subset.pdf)
[//]: # (Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation.R \)
[//]: # (  $path_grounded_pq_touching_reliable_tsv_gz \)
[//]: # (  0 25000000 \)
[//]: # (  91 0.8 \)
[//]: # (  0 \)
[//]: # (  1 1 \)
[//]: # (  $i \)
[//]: # (  0.9 \)
[//]: # (  <&#40;zgrep '^chm\|grch\|^HG002#' $path_grounded_pq_touching_reliable_tsv_gz | cut -f 1 | grep "chr$i\|^HG002" | sort | uniq&#41; \)
[//]: # (  /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \)
[//]: # (  /lizardfs/guarracino/chromosome_communities/untangle/grounded/$PREFIX.n1.HG002.pdf)
[//]: # ()
[//]: # (i=14)
[//]: # (path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chr${i}.e$e.m$m.grounded.pq_touching.reliable.tsv.gz)
[//]: # (PREFIX=$&#40;basename $path_grounded_pq_touching_reliable_tsv_gz .tsv.gz&#41;;)
[//]: # (Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation.R \)
[//]: # (  $path_grounded_pq_touching_reliable_tsv_gz \)
[//]: # (  0 25000000 \)
[//]: # (  91 0.8 \)
[//]: # (  0 \)
[//]: # (  1 1 \)
[//]: # (  $i \)
[//]: # (  0.9 \)
[//]: # (  <&#40;echo chm13#chr14 grch38#chr14 HG002#MAT#chr14.prox HG002#PAT#chr14.prox HG00735#1#JAHBCH010000039.1 HG00741#2#JAHALX010000038.1 HG01978#1#JAGYVS010000055.1 HG02630#1#JAHAOQ010000067.1 | tr ' ' '\n'&#41; \)
[//]: # (  /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \)
[//]: # (  ~/$PREFIX.n1.subset.pdf)
[//]: # (Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation.R \)
[//]: # (  $path_grounded_pq_touching_reliable_tsv_gz \)
[//]: # (  0 25000000 \)
[//]: # (  91 0.8 \)
[//]: # (  0.8 \)
[//]: # (  5 1 \)
[//]: # (  $i \)
[//]: # (  0.9 \)
[//]: # (  <&#40;echo chm13#chr14 grch38#chr14 HG002#MAT#chr14.prox HG002#PAT#chr14.prox HG00735#1#JAHBCH010000039.1 HG00741#2#JAHALX010000038.1 HG01978#1#JAGYVS010000055.1 HG02630#1#JAHAOQ010000067.1 | tr ' ' '\n'&#41; \)
[//]: # (  /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \)
[//]: # (  /lizardfs/guarracino/chromosome_communities/untangle/grounded/$PREFIX.n5.subset.pdf)
[//]: # (Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation.R \)
[//]: # (  $path_grounded_pq_touching_reliable_tsv_gz \)
[//]: # (  0 25000000 \)
[//]: # (  91 0.8 \)
[//]: # (  0 \)
[//]: # (  1 1 \)
[//]: # (  $i \)
[//]: # (  0.9 \)
[//]: # (  <&#40;zgrep '^chm\|grch\|^HG002#' $path_grounded_pq_touching_reliable_tsv_gz | cut -f 1 | grep "chr$i\|^HG002" | sort | uniq&#41; \)
[//]: # (  /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \)
[//]: # (  /lizardfs/guarracino/chromosome_communities/untangle/grounded/$PREFIX.n1.HG002.pdf)
[//]: # ()
[//]: # (i=15)
[//]: # (path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chr${i}.e$e.m$m.grounded.pq_touching.reliable.tsv.gz)
[//]: # (PREFIX=$&#40;basename $path_grounded_pq_touching_reliable_tsv_gz .tsv.gz&#41;;)
[//]: # (Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation.R \)
[//]: # (  $path_grounded_pq_touching_reliable_tsv_gz \)
[//]: # (  0 25000000 \)
[//]: # (  91 0.8 \)
[//]: # (  0 \)
[//]: # (  1 1 \)
[//]: # (  $i \)
[//]: # (  0.9 \)
[//]: # (  <&#40;echo chm13#chr15 grch38#chr15 HG002#MAT#chr15.prox HG002#PAT#chr15.prox HG00741#2#JAHALX010000004.1 HG02486#2#JAGYVL010000058.1 HG03486#2#JAHEOP010000088.1 NA18906#2#JAHEON010000012.1 | tr ' ' '\n'&#41; \)
[//]: # (  /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \)
[//]: # (  ~/$PREFIX.n1.subset.pdf)
[//]: # (Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation.R \)
[//]: # (  $path_grounded_pq_touching_reliable_tsv_gz \)
[//]: # (  0 25000000 \)
[//]: # (  91 0.8 \)
[//]: # (  0.8 \)
[//]: # (  5 1 \)
[//]: # (  $i \)
[//]: # (  0.9 \)
[//]: # (  <&#40;echo chm13#chr15 grch38#chr15 HG002#MAT#chr15.prox HG002#PAT#chr15.prox HG00741#2#JAHALX010000004.1 HG02486#2#JAGYVL010000058.1 HG03486#2#JAHEOP010000088.1 NA18906#2#JAHEON010000012.1 | tr ' ' '\n'&#41; \)
[//]: # (  /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \)
[//]: # (  /lizardfs/guarracino/chromosome_communities/untangle/grounded/$PREFIX.n5.subset.pdf)
[//]: # (Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation.R \)
[//]: # (  $path_grounded_pq_touching_reliable_tsv_gz \)
[//]: # (  0 25000000 \)
[//]: # (  91 0.8 \)
[//]: # (  0 \)
[//]: # (  1 1 \)
[//]: # (  $i \)
[//]: # (  0.9 \)
[//]: # (  <&#40;zgrep '^chm\|grch\|^HG002#' $path_grounded_pq_touching_reliable_tsv_gz | cut -f 1 | grep "chr$i\|^HG002" | sort | uniq&#41; \)
[//]: # (  /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \)
[//]: # (  /lizardfs/guarracino/chromosome_communities/untangle/grounded/$PREFIX.n1.HG002.pdf)
[//]: # ()
[//]: # (i=21)
[//]: # (path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chr${i}.e$e.m$m.grounded.pq_touching.reliable.tsv.gz)
[//]: # (PREFIX=$&#40;basename $path_grounded_pq_touching_reliable_tsv_gz .tsv.gz&#41;;)
[//]: # (Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation.R \)
[//]: # (  $path_grounded_pq_touching_reliable_tsv_gz \)
[//]: # (  0 25000000 \)
[//]: # (  91 0.8 \)
[//]: # (  0 \)
[//]: # (  1 1 \)
[//]: # (  $i \)
[//]: # (  0.9 \)
[//]: # (  <&#40;echo chm13#chr21 grch38#chr21 HG002#MAT#chr21.prox HG002#PAT#chr21.prox HG00735#2#JAHBCG010000066.1 HG02886#1#JAHAOU010000106.1 NA18906#1#JAHEOO010000072.1 NA19240#2#JAHEOL010000065.1 | tr ' ' '\n'&#41; \)
[//]: # (  /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \)
[//]: # (  ~/$PREFIX.n1.subset.pdf)
[//]: # (Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation.R \)
[//]: # (  $path_grounded_pq_touching_reliable_tsv_gz \)
[//]: # (  0 25000000 \)
[//]: # (  91 0.8 \)
[//]: # (  0 \)
[//]: # (  1 1 \)
[//]: # (  $i \)
[//]: # (  0.9 \)
[//]: # (  <&#40;zgrep '^chm\|grch\|^HG002#' $path_grounded_pq_touching_reliable_tsv_gz | cut -f 1 | grep "chr$i\|^HG002" | sort | uniq&#41; \)
[//]: # (  /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \)
[//]: # (  /lizardfs/guarracino/chromosome_communities/untangle/grounded/$PREFIX.n1.HG002.pdf)
[//]: # (  )
[//]: # (i=22)
[//]: # (path_grounded_pq_touching_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chr${i}.e$e.m$m.grounded.pq_touching.reliable.tsv.gz)
[//]: # (PREFIX=$&#40;basename $path_grounded_pq_touching_reliable_tsv_gz .tsv.gz&#41;;)
[//]: # (Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation.R \)
[//]: # (  $path_grounded_pq_touching_reliable_tsv_gz \)
[//]: # (  0 25000000 \)
[//]: # (  91 0.8 \)
[//]: # (  0 \)
[//]: # (  1 1 \)
[//]: # (  $i \)
[//]: # (  0.9 \)
[//]: # (  <&#40;echo chm13#chr22 grch38#chr22 HG002#MAT#chr22.prox HG002#PAT#chr22.prox HG00735#1#JAHBCH010000040.1 HG01361#1#JAGYYX010000045.1 HG02055#1#JAHEPK010000087.1 HG03098#1#JAHEPM010000147.1 | tr ' ' '\n'&#41; \)
[//]: # (  /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \)
[//]: # (  ~/$PREFIX.n1.subset.pdf)
[//]: # (Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_with_annotation.R \)
[//]: # (  $path_grounded_pq_touching_reliable_tsv_gz \)
[//]: # (  0 25000000 \)
[//]: # (  91 0.8 \)
[//]: # (  0 \)
[//]: # (  1 1 \)
[//]: # (  $i \)
[//]: # (  0.9 \)
[//]: # (  <&#40;zgrep '^chm\|grch\|^HG002#' $path_grounded_pq_touching_reliable_tsv_gz | cut -f 1 | grep "chr$i\|^HG002" | sort | uniq&#41; \)
[//]: # (  /lizardfs/guarracino/chromosome_communities/data/annotation/hgt_genome_euro_chr${i}_0_25Mbp.png \)
[//]: # (  /lizardfs/guarracino/chromosome_communities/untangle/grounded/$PREFIX.n1.HG002.pdf)
[//]: # (```)
[//]: # (Merged plots &#40;not used&#41;:)
[//]: # ()
[//]: # (```shell)
[//]: # (#for e in 50000 ; do)
[//]: # (#  for m in 1000 ; do)
[//]: # (#    echo "-e $e -m $m")
[//]: # (#)
[//]: # (#    path_grounded_pq_touching_reliable_all_chromosomes_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.chm13#chrACRO.e$e.m$m.grounded.pq_touching.reliable.tsv.gz)
[//]: # (#)
[//]: # (#    Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_all2.R ${path_grounded_pq_touching_reliable_all_chromosomes_tsv_gz} "-e $e -m $m" 0 25000000 120 200 1 1)
[//]: # (#  done)
[//]: # (#done)
[//]: # (#)
[//]: # (#mv /lizardfs/guarracino/chromosome_communities/untangle/grounded/*.pdf /lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/)
[//]: # ()
[//]: # (# Merge)
[//]: # (#for e in 50000 ; do)
[//]: # (#  for m in 1000 ; do)
[//]: # (#    for j in 0.8 0.95; do)
[//]: # (#      echo "-e $e -m $m -j $j")
[//]: # (#        )
[//]: # (#      j_str=$&#40;echo $j | sed 's/\.//g'&#41;)
[//]: # (#      PDFs=$&#40;&#40;seq 1 5; seq 10 10 50&#41; | while read n; do \)
[//]: # (#        echo /lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/$prefix.untangle.chm13#chrACRO.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.reliable.tsv.gz.pdf)
[//]: # (#      done | tr '\n' ' '&#41;)
[//]: # (#      #echo $PDFs)
[//]: # (#        )
[//]: # (#      # Too slow)
[//]: # (#      #gs -sDEVICE=pdfwrite \)
[//]: # (#      #  -dNOPAUSE -dBATCH -dSAFER \)
[//]: # (#      #  -sOutputFile=/lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/$prefix.untangle.chm13#chrACRO.e$e.m$m.j${j_str}.grounded.pq_touching.tsv.gz.pdf \)
[//]: # (#      #  $PDFs)
[//]: # (#        )
[//]: # (#      /gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite $PDFs /lizardfs/guarracino/chromosome_communities/untangle/grounded/pdf/$prefix.untangle.chm13#chrACRO.e$e.m$m.j${j_str}.merged.grounded.pq_touching.reliable.tsv.gz.pdf)
[//]: # (#    done)
[//]: # (#  done)
[//]: # (#done)
[//]: # (```)
[//]: # ()
[//]: # (Single plots &#40;not used&#41;:)
[//]: # ()
[//]: # (```shell)
[//]: # (#for e in 5000 50000 100000; do)
[//]: # (#  for m in 500 1000 10000; do)
[//]: # (#      for j in 0 0.8; do)
[//]: # (#        j_str=$&#40;echo $j | sed 's/\.//g'&#41;)
[//]: # (#        &#40;seq 1 5; seq 10 10 50&#41; | while read n; do )
[//]: # (#            echo "-e $e -m $m -j $j -n $n")
[//]: # (#    )
[//]: # (#            grep chm13 /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.1kbps.pq_contigs.union.hg002prox.fa.gz.fai | cut -f 1 | while read ref; do)
[//]: # (#              echo $ref)
[//]: # (#              )
[//]: # (#              path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz)
[//]: # (#            )
[//]: # (#              Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle.R ${path_grounded_pq_touching_tsv_gz} "$ref -e $e -m $m -j $j -n $n")
[//]: # (#            done)
[//]: # (#        done)
[//]: # (#      done)
[//]: # (#    done)
[//]: # (#done)
[//]: # (```)
[//]: # ()
[//]: # ()
[//]: # (#### Plot Mobin's annotations &#40;not used&#41;)
[//]: # ()
[//]: # (Take only reliable blocks &#40;flagged with "Hh" or "Hc", https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/README.md#components&#41;)
[//]: # ()
[//]: # (```shell)
[//]: # (e=50000)
[//]: # (m=1000)
[//]: # (j=0.8)
[//]: # (j_str=$&#40;echo $j | sed 's/\.//g'&#41;)
[//]: # (n=1)
[//]: # (ref=chm13#chr13)
[//]: # ()
[//]: # (# Brutal: remove untangled regions if they overlap with unreliable regions)
[//]: # (touch xyz.tsv)
[//]: # (cat $path_targets_txt | while read ref; do)
[//]: # (    echo $ref)
[//]: # (    path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz)
[//]: # (    )
[//]: # (    touch z.tsv)
[//]: # (    zcat $path_grounded_pq_touching_tsv_gz | sed '1d' | grep 'HG002#MAT\|HG002#PAT\|HG01978#MAT\|HG01978#PAT\|bakeoff' -v | cut -f 1 | sort | uniq | while read CONTIG; do)
[//]: # (      SAMPLE=$&#40; echo $CONTIG | cut -f 1 -d '#'&#41;)
[//]: # (    )
[//]: # (      path_unreliable_bed=/lizardfs/erikg/HPRC/year1v2genbank/annotations/unreliable/$SAMPLE.hifi.flagger_final.bed)
[//]: # (      if [[ -s $path_unreliable_bed ]]; then)
[//]: # (        echo $CONTIG "--->" $SAMPLE)
[//]: # (        )
[//]: # (        zgrep "^$CONTIG" $path_grounded_pq_touching_tsv_gz | sort -k 2n | awk -v OFS='\t' '{print&#40;$1,$2,$3,$4"_"$5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11"_"$12"_"$13&#41;}' > x.bed)
[//]: # (        grep $CONTIG $path_unreliable_bed > y.bed)
[//]: # (        # -A: remove entire feature if any overlap)
[//]: # (        bedtools subtract -a x.bed -b y.bed -A |\)
[//]: # (          awk -v OFS='\t' '{split&#40;$4, a, "_"&#41;; print&#40;$1,$2,$3,a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10]&#41;}' >> xyz.tsv)
[//]: # (      else)
[//]: # (        zgrep "^$CONTIG" $path_grounded_pq_touching_tsv_gz | sort -k 2n >> xyz.tsv)
[//]: # (      fi)
[//]: # (    done)
[//]: # (    )
[//]: # (    zcat $path_grounded_pq_touching_tsv_gz | sed '1d' | grep 'HG002#MAT\|HG002#PAT\|HG01978#MAT\|HG01978#PAT\|bakeoff' >> xyz.tsv)
[//]: # (done)
[//]: # ()
[//]: # (path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz)
[//]: # (cat \)
[//]: # (  <&#40; zcat $path_grounded_pq_touching_tsv_gz | head -n 1 &#41; \)
[//]: # (   xyz.tsv | pigz -c > xyz.tsv.gz)
[//]: # (rm xyz.tsv)
[//]: # ()
[//]: # (Rscript scripts/plot_untangle_all.R xyz.tsv.gz "-e 50000 -m 1000 -j 0.8 -n 1" 0 25000000 360 200)
[//]: # ()
[//]: # ()
[//]: # (# Plot additional tracks for the annotations)
[//]: # (touch xyz.tsv)
[//]: # (cat $path_targets_txt | while read ref; do)
[//]: # (    echo $ref)
[//]: # (    path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz)
[//]: # (    )
[//]: # (    touch z.tsv)
[//]: # (    zcat $path_grounded_pq_touching_tsv_gz | sed '1d' | grep 'HG002#MAT\|HG002#PAT\|HG01978#MAT\|HG01978#PAT\|bakeoff' -v | cut -f 1 | sort | uniq | while read CONTIG; do)
[//]: # (      SAMPLE=$&#40; echo $CONTIG | cut -f 1 -d '#'&#41;)
[//]: # (    )
[//]: # (      path_unreliable_bed=/lizardfs/erikg/HPRC/year1v2genbank/annotations/unreliable/$SAMPLE.hifi.flagger_final.bed)
[//]: # (      if [[ -s $path_unreliable_bed ]]; then)
[//]: # (        echo $CONTIG "--->" $SAMPLE)
[//]: # (        )
[//]: # (        zgrep "^$CONTIG" $path_grounded_pq_touching_tsv_gz | sort -k 2n | awk -v OFS='\t' '{print&#40;$1,$2,$3,$4"_"$5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11"_"$12"_"$13&#41;}' > x.bed)
[//]: # (        grep $CONTIG $path_unreliable_bed > y.bed)
[//]: # (        )
[//]: # (        # wao: write the original A and B entries plus the number of base pairs of overlap between the two features.)
[//]: # (        bedtools intersect -a x.bed -b y.bed -wo >> z.tsv)
[//]: # (      fi)
[//]: # (    done)
[//]: # (    )
[//]: # (    cat \)
[//]: # (      <&#40; zcat $path_grounded_pq_touching_tsv_gz | sed '1d' &#41; \)
[//]: # (      <&#40; python3 scripts/get_annotation_track.py z.tsv &#41; | tr ' ' '\t' >> xyz.tsv)
[//]: # ()
[//]: # (    rm z.tsv)
[//]: # (done)
[//]: # ()
[//]: # (path_grounded_pq_touching_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle/grounded/$prefix.untangle.$ref.e$e.m$m.j${j_str}.n$n.grounded.pq_touching.tsv.gz)
[//]: # (cat \)
[//]: # (  <&#40; zcat $path_grounded_pq_touching_tsv_gz | head -n 1 &#41; \)
[//]: # (   xyz.tsv | pigz -c > xyz.annot.tsv.gz)
[//]: # (rm xyz.tsv)
[//]: # ()
[//]: # (Rscript scripts/plot_untangle_all.R xyz.annot.tsv.gz "-e 50000 -m 1000 -j 0.8 -n 1" 0 25000000 120 200)
[//]: # (```)


## Variant calling

Call variants in a haploid setting:

```shell
path_input_gfa=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.gfa
path_chm13_vcf_gz=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.chm13.haploid.vcf.gz
sbatch -p headnode -c 48 --job-name vgchm13 --wrap "\time -v vg deconstruct -P chm13 -H '?' -e -a --ploidy 1 -t 48 $path_input_gfa | bgzip -@ 48 -c > $path_chm13_vcf_gz && tabix $path_chm13_vcf_gz"

# Take SNVs
zcat chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.chm13.haploid.vcf.gz | awk -F '\t' '($0 ~ /^#/ || (length($4)==1 && length($5)==1))' | bgzip -c -@ 48 > chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.chm13.haploid.snv.vcf.gz

# Take SNPs and normalize
f=chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.chm13.haploid
vcftools --gzvcf $f.vcf.gz --remove-indels --recode --recode-INFO-all --stdout |\
  bcftools norm -f /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz -c s -m - \
  > $f.snv.norm.vcf
# Remove HG002-HiFi contigs
vcfsamplenames $f.snv.norm.vcf | grep -v HG002#1 | grep -v HG002#2 > keep.samples
vcfkeepsamples $f.snv.norm.vcf $(cat keep.samples) | vcffixup - | bgzip -c -@ 48 > $f.snv.norm.no_HG002-hifi.vcf.gz


## acro-contigs
#path_input_gfa=/lizardfs/guarracino/chromosome_communities/graphs/chrA.pan+HG002chrAprox.s100k.l300k.p98.n100/chrA.pan+HG002chrAprox.fa.gz.952ce3b.4030258.0e65f78.smooth.fix.gfa
#path_input_sed_gfa=/lizardfs/guarracino/chromosome_communities/graphs/chrA.pan+HG002chrAprox.s100k.l300k.p98.n100/chrA.pan+HG002chrAprox.fa.gz.952ce3b.4030258.0e65f78.smooth.fix.sed.gfa
#path_chm13_vcf_gz=/lizardfs/guarracino/chromosome_communities/graphs/chrA.pan+HG002chrAprox.s100k.l300k.p98.n100/chrA.pan+HG002chrAprox.fa.gz.952ce3b.4030258.0e65f78.smooth.fix.sed.c1000.chm13.vcf.gz
#
## sed replaces only the first instance on a line by default (without the /g modifier)
## To have names like NA21309-1#1#JAHEPC010000450.1 and call haploid genotypes with -H
#sed 's/#/-/' $path_input_gfa | sed 's/#/#1#/' > $path_input_sed_gfa
#sbatch -p headnode -c 48 --job-name vgchm13 --wrap '\time -v vg deconstruct -P chm13 -H "#" -e -a -c 1000 -t 48 '$path_input_sed_gfa' | bgzip -@ 48 -c > '$path_chm13_vcf_gz' && tabix '$path_chm13_vcf_gz
```
