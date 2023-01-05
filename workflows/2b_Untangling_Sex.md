# Untangling

## Sex chromosome

Use HG002's chrY as an alternative reference, as GRCh38's chrY is incomplete. Include also HG002's chrX.

```shell
cd /lizardfs/guarracino/chromosome_communities/assemblies

# Get HG002 chromosome X and Y

PATH_CHM13_FA_GZ=/lizardfs/erikg/human/chm13v2.0.fa.gz
PATH_GRCH38_FA_GZ=/lizardfs/guarracino/chromosome_communities/assemblies/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2_maskedcentromeres_auto_XY.fa.gz

PATH_HPRCY1_FA_GZ=/lizardfs/guarracino/chromosome_communities/assemblies/HPRCy1v2genbank.fa.gz

cat /lizardfs/guarracino/chromosome_communities/assemblies/partitioning/*.partitions.tsv | awk -v OFS='\t' '($2 == "chrX" || $2 == "chrY") {print $1}' | sort > chrXY_contigs.txt

# HiFi-based assemblies + hg01978 verkko + hg002-bakeoff + hg002 verkko
cat \
  <( samtools faidx $PATH_CHM13_FA_GZ $(grep 'chrX\|chrY' $PATH_CHM13_FA_GZ.fai | cut -f 1) | sed 's/^>/>chm13#/g' ) \
  <( samtools faidx $PATH_GRCH38_FA_GZ $(grep 'chrX\|chrY' $PATH_GRCH38_FA_GZ.fai | grep '_' -v | cut -f 1) ) \
  <(samtools faidx $PATH_HPRCY1_FA_GZ $(comm -12 <(cut -f 1 "$PATH_HPRCY1_FA_GZ".fai | sort) <(cat chrXY_contigs.txt))) \
  <(samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/hg01978.mat.fa $(comm -12 <(cut -f 1 /lizardfs/guarracino/chromosome_communities/assemblies/hg01978.mat.fa.fai | sort) <(cat chrXY_contigs.txt))) \
  <(samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/hg01978.pat.fa $(comm -12 <(cut -f 1 /lizardfs/guarracino/chromosome_communities/assemblies/hg01978.pat.fa.fai | sort) <(cat chrXY_contigs.txt))) \
  <(samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/hg002-bakeoff.mat.fa $(comm -12 <(cut -f 1 /lizardfs/guarracino/chromosome_communities/assemblies/hg002-bakeoff.mat.fa.fai | sort) <(cat chrXY_contigs.txt))) \
  <(samtools faidx /lizardfs/guarracino/chromosome_communities/assemblies/hg002-bakeoff.pat.fa $(comm -12 <(cut -f 1 /lizardfs/guarracino/chromosome_communities/assemblies/hg002-bakeoff.pat.fa.fai | sort) <(cat chrXY_contigs.txt))) | \
  bgzip -@ 48 -c > chrSEX+refs.fa.gz
samtools faidx chrSEX+refs.fa.gz
```


### Pangenome building

Apply `pggb` on the chromosome-partitioned HPRC dataset:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/graphs

RUN_PGGB=/home/guarracino/tools/pggb/pggb-a4a6668d9ece42c80ce69dc354f0cb59a849286f

num_of_haplotypes=$(cut -f 1,2 -d '#' /lizardfs/guarracino/chromosome_communities/assemblies/chrSEX+refs.fa.gz.fai | sort | uniq | wc -l)

sbatch -p highmem -c 48 --job-name sexpggb --wrap "hostname; cd /scratch && $RUN_PGGB -i /lizardfs/guarracino/chromosome_communities/assemblies/chrSEX+refs.fa.gz -o chrSEX+refs.s50k.l250k.p98.n${num_of_haplotypes} -t 48 -s 50k -l 250k -p 98 -n ${num_of_haplotypes} -k 311 -G 13117,13219 -O 0.03 -T 48 -v -V chm13:#,grch38:#; mv /scratch/chrSEX+refs.s50k.l250k.p98.n${num_of_haplotypes} /lizardfs/guarracino/chromosome_communities/graphs";

sbatch -p highmem -c 48 --job-name sexpggb --wrap "hostname; cd /scratch && $RUN_PGGB -i /lizardfs/guarracino/chromosome_communities/assemblies/chrSEX+refs.fa.gz -o chrSEX+refs.s50k.l250k.p90.n${num_of_haplotypes} -t 48 -s 50k -l 250k -p 90 -n ${num_of_haplotypes} -k 311 -G 13117,13219 -O 0.03 -T 48 -v -V chm13:#,grch38:#; mv /scratch/chrSEX+refs.s50k.l250k.p90.n${num_of_haplotypes} /lizardfs/guarracino/chromosome_communities/graphs";
```


### Untangling

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/untangle_sex/

path_targets_txt=/lizardfs/guarracino/chromosome_communities/untangle_sex/chm13.SEX.target_paths.txt
path_fasta_fai=/lizardfs/guarracino/chromosome_communities/assemblies/chrSEX+refs.fa.gz.fai
grep chm13 $path_fasta_fai | cut -f 1 > $path_targets_txt


path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrSEX+refs.s50k.l250k.p98.n102/chrSEX+refs.fa.gz.2ed2c67.04f1c29.22fc5c8.smooth.final.og
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrSEX+refs.s50k.l250k.p90.n102/chrSEX+refs.fa.gz.d1ae18e.04f1c29.d09fe3b.smooth.final.og
prefix=$(basename $path_input_og .og)

RUN_ODGI=/home/guarracino/tools/odgi/bin/odgi-454197fa29b772050c3135d5de47c816ce38e62c
```


Untangle with respect to all sex chromosomes and emit the cut points:

```shell
for e in 50000; do
  for m in 1000; do
    path_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle_sex/$prefix.untangle.chm13#SEX.e$e.m$m.j0.n100.bed.gz
    path_cut_points_txt=/lizardfs/guarracino/chromosome_communities/untangle_sex/$prefix.untangle.chm13#SEX.e$e.m$m.j0.n100.cut_points.txt
    
    if [[ ! -s ${path_cut_points_txt} ]]; then
      echo "-e $e -m $m"
      sbatch -p workers -c 48 --job-name sexuntangle --wrap "\time -v $RUN_ODGI untangle -t 48 -P -i $path_input_og -R $path_targets_txt -e $e -m $m --cut-points-output $path_cut_points_txt -j 0 -n 100 | pigz -c > $path_bed_gz"
    fi;
  done
done
```


Fix best hits (if there are multiple best hits, put as first the target-chromosome of origin of the contig):

```shell
cat \
            <(cat assemblies/partitioning/*.partitions.tsv | sed 's/chr/chm13#chr/') \
            <(echo -e HG002-bakeoff#PAT#SY_unloc_13"\t"chm13#chrY) \
            <(echo -e HG002-bakeoff#PAT#SY_unloc_15"\t"chm13#chrY) \
            <(echo -e HG002-bakeoff#PAT#SY_unloc_19"\t"chm13#chrY) \
            <(echo -e HG002-bakeoff#PAT#SY_unloc_21"\t"chm13#chrY) \
            <(echo -e HG002-bakeoff#PAT#SY_unloc_22"\t"chm13#chrY) \
            <(echo -e HG002-bakeoff#PAT#SY_unloc_22"\t"chm13#chrY) \
            <(echo -e HG002#1#JAHKSE010000383.1"\t"chm13#chrY) \
            <(echo -e HG002#1#JAHKSE010000141.1"\t"chm13#chrY) \
            <(echo -e HG002#1#JAHKSE010000378.1"\t"chm13#chrY) \
            <(echo -e HG002#1#JAHKSE010000376.1"\t"chm13#chrY) \
            <(echo -e HG002#1#JAHKSE010000341.1"\t"chm13#chrY) \
            <(echo -e HG002#1#JAHKSE010000447.1"\t"chm13#chrY) \
            <(echo -e HG005#1#JAHEPO010000097.1"\t"chm13#chrY) \
            <(echo -e HG005#1#JAHEPO010000270.1"\t"chm13#chrY) \
            <(echo -e HG00621#1#JAHBCD010000210.1"\t"chm13#chrY) \
            <(echo -e HG005#1#JAHEPO010000116.1"\t"chm13#chrY) \
            <(echo -e HG005#1#JAHEPO010000199.1"\t"chm13#chrX) \
            <(echo -e HG00673#2#JAHBBY010000245.1"\t"chm13#chrX) \
            <(echo -e HG005#2#JAHEPN010000186.1"\t"chm13#chrX) \
            <(echo -e HG00733#1#JAHEPQ010000312.1"\t"chm13#chrX) \
            <(echo -e HG01071#2#JAHBCE010000198.1"\t"chm13#chrX) \
            <(echo -e HG01071#1#JAHBCF010000145.1"\t"chm13#chrX) \
            <(echo -e HG01109#1#JAHEPA010000188.1"\t"chm13#chrY) \
            <(echo -e HG01109#1#JAHEPA010000268.1"\t"chm13#chrY) \
            <(echo -e HG01109#1#JAHEPA010000227.1"\t"chm13#chrY) \
            <(echo -e HG01109#1#JAHEPA010000307.1"\t"chm13#chrY) \
            <(echo -e HG01109#1#JAHEPA010000343.1"\t"chm13#chrY) \
            <(echo -e HG01109#1#JAHEPA010000333.1"\t"chm13#chrY) \
            <(echo -e HG01109#1#JAHEPA010000405.1"\t"chm13#chrY) \
            <(echo -e HG01109#1#JAHEPA010000337.1"\t"chm13#chrY) \
            <(echo -e HG01109#1#JAHEPA010000344.1"\t"chm13#chrY) \
            <(echo -e HG01109#1#JAHEPA010000442.1"\t"chm13#chrY) \
            <(echo -e HG01109#1#JAHEPA010000458.1"\t"chm13#chrY) \
            <(echo -e HG01109#1#JAHEPA010000325.1"\t"chm13#chrY) \
            <(echo -e HG01123#2#JAGYYY010000193.1"\t"chm13#chrX) \
            <(echo -e HG01175#1#JAHAMA010000276.1"\t"chm13#chrX) \
            <(echo -e HG01175#2#JAHALZ010000239.1"\t"chm13#chrX) \
            <(echo -e HG01243#1#JAHEOY010000277.1"\t"chm13#chrY) \
            <(echo -e HG01243#1#JAHEOY010000307.1"\t"chm13#chrY) \
            <(echo -e HG01243#1#JAHEOY010000291.1"\t"chm13#chrY) \
            <(echo -e HG01258#1#JAGYYV010000281.1"\t"chm13#chrY) \
            <(echo -e HG01358#2#JAGYZA010000069.1"\t"chm13#chrX) \
            <(echo -e HG01891#2#JAGYVN010000149.1"\t"chm13#chrX) \
            <(echo -e HG01891#1#JAGYVO010000131.1"\t"chm13#chrX) \
            <(echo -e HG01928#1#JAGYVQ010000246.1"\t"chm13#chrY) \
            <(echo -e HG01952#1#JAHAME010000217.1"\t"chm13#chrY) \
            <(echo -e HG02055#1#JAHEPK010000284.1"\t"chm13#chrY) \
            <(echo -e HG02145#1#JAHKSG010000213.1"\t"chm13#chrY) \
            <(echo -e HG02145#1#JAHKSG010000422.1"\t"chm13#chrY) \
            <(echo -e HG02145#1#JAHKSG010000527.1"\t"chm13#chrY) \
            <(echo -e HG02145#1#JAHKSG010000498.1"\t"chm13#chrY) \
            <(echo -e HG02145#1#JAHKSG010000544.1"\t"chm13#chrY) \
            <(echo -e HG02145#1#JAHKSG010000728.1"\t"chm13#chrY) \
            <(echo -e HG02145#1#JAHKSG010000766.1"\t"chm13#chrY) \
            <(echo -e HG02145#1#JAHKSG010000668.1"\t"chm13#chrY) \
            <(echo -e HG02148#1#JAHAMG010000114.1"\t"chm13#chrX) \
            <(echo -e HG02486#1#JAGYVM010000047.1"\t"chm13#chrY) \
            <(echo -e HG02148#2#JAHAMF010000125.1"\t"chm13#chrX) \
            <(echo -e HG02486#1#JAGYVM010000189.1"\t"chm13#chrY) \
            <(echo -e HG02572#1#JAHAOW010000422.1"\t"chm13#chrY) \
            <(echo -e HG02572#1#JAHAOW010000309.1"\t"chm13#chrY) \
            <(echo -e HG02572#1#JAHAOW010000471.1"\t"chm13#chrY) \
            <(echo -e HG02572#1#JAHAOW010000523.1"\t"chm13#chrY) \
            <(echo -e HG02572#1#JAHAOW010000458.1"\t"chm13#chrY) \
            <(echo -e HG02572#1#JAHAOW010000535.1"\t"chm13#chrY) \
            <(echo -e HG02572#2#JAHAOV010000449.1"\t"chm13#chrY) \
            <(echo -e HG02622#2#JAHAON010000142.1"\t"chm13#chrX) \
            <(echo -e HG02717#1#JAHAOS010000199.1"\t"chm13#chrY) \
            <(echo -e HG02886#2#JAHAOT010000422.1"\t"chm13#chrX) \
            <(echo -e HG03098#1#JAHEPM010000129.1"\t"chm13#chrY) \
            <(echo -e HG03098#1#JAHEPM010000235.1"\t"chm13#chrY) \
            <(echo -e HG03098#1#JAHEPM010000398.1"\t"chm13#chrY) \
            <(echo -e HG03098#1#JAHEPM010000264.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000175.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000301.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000359.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000353.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000402.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000380.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000404.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000442.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000397.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000407.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000414.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000466.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000510.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000522.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000516.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000418.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000548.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000476.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000591.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000542.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000554.1"\t"chm13#chrY) \
            <(echo -e HG03516#2#JAGYYS010000236.1"\t"chm13#chrX) \
            <(echo -e HG03579#1#JAGYVU010000286.1"\t"chm13#chrY) \
            <(echo -e HG03579#1#JAGYVU010000293.1"\t"chm13#chrY) \
            <(echo -e HG03579#1#JAGYVU010000236.1"\t"chm13#chrY) \
            <(echo -e HG03579#1#JAGYVU010000411.1"\t"chm13#chrY) \
            <(echo -e HG03579#1#JAGYVU010000348.1"\t"chm13#chrY) \
            <(echo -e HG03579#1#JAGYVU010000340.1"\t"chm13#chrY) \
            <(echo -e HG03579#1#JAGYVU010000420.1"\t"chm13#chrY) \
            <(echo -e HG03540#2#JAGYVX010000220.1"\t"chm13#chrX) \
            <(echo -e NA18906#2#JAHEON010000121.1"\t"chm13#chrX) \
            <(echo -e NA19240#1#JAHEOM010000255.1"\t"chm13#chrX) \
            <(echo -e NA21309#1#JAHEPC010000444.1"\t"chm13#chrX) \
            <(echo -e NA20129#2#JAHEPD010000373.1"\t"chm13#chrX) \
            <(echo -e NA21309#1#JAHEPC010000484.1"\t"chm13#chrX) \
            <(echo -e NA21309#1#JAHEPC010000136.1"\t"chm13#chrX) \
            <(echo -e HG005#1#JAHEPO010000137.1"\t"chm13#chrY) \
            <(echo -e HG005#1#JAHEPO010000240.1"\t"chm13#chrY) \
            <(echo -e HG005#1#JAHEPO010000251.1"\t"chm13#chrY) \
            <(echo -e HG005#1#JAHEPO010000358.1"\t"chm13#chrY) \
            <(echo -e HG005#1#JAHEPO010000287.1"\t"chm13#chrX) \
            <(echo -e HG005#2#JAHEPN010000184.1"\t"chm13#chrX) \
            <(echo -e HG005#2#JAHEPN010000184.1"\t"chm13#chrX) \
            <(echo -e HG005#1#JAHEPO010000150.1"\t"chm13#chrY) \
            <(echo -e HG002#1#JAHKSE010000138.1"\t"chm13#chrY) \
            <(echo -e HG005#1#JAHEPO010000129.1"\t"chm13#chrY) \
            <(echo -e HG005#1#JAHEPO010000146.1"\t"chm13#chrY) \
            <(echo -e HG00673#1#JAHBBZ010000301.1"\t"chm13#chrY) \
            <(echo -e HG00673#1#JAHBBZ010000306.1"\t"chm13#chrY) \
            <(echo -e HG00673#2#JAHBBY010000140.1"\t"chm13#chrX) \
            <(echo -e HG00673#2#JAHBBY010000140.1"\t"chm13#chrX) \
            <(echo -e HHG00733#1#JAHEPQ010000190.1"\t"chm13#chrX) \
            <(echo -e HG00673#2#JAHBBY010000210.1"\t"chm13#chrX) \
            <(echo -e HG00733#1#JAHEPQ010000190.1"\t"chm13#chrX) \
            <(echo -e HG00733#2#JAHEPP010000392.1"\t"chm13#chrX) \
            <(echo -e HG005#1#JAHEPO010000178.1"\t"chm13#chrY) \
            <(echo -e HG00673#1#JAHBBZ010000257.1"\t"chm13#chrY) \
            <(echo -e HG01109#1#JAHEPA010000218.1"\t"chm13#chrY) \
            <(echo -e HG01109#1#JAHEPA010000263.1"\t"chm13#chrY) \
            <(echo -e HG01109#1#JAHEPA010000278.1"\t"chm13#chrY) \
            <(echo -e HG01109#1#JAHEPA010000428.1"\t"chm13#chrY) \
            <(echo -e G01123#1#JAGYYZ010000259.1"\t"chm13#chrX) \
            <(echo -e HG01175#2#JAHALZ010000133.1"\t"chm13#chrX) \
            <(echo -e HG01243#1#JAHEOY010000179.1"\t"chm13#chrY) \
            <(echo -e HG01258#1#JAGYYV010000194.1"\t"chm13#chrY) \
            <(echo -e HG01952#1#JAHAME010000216.1"\t"chm13#chrY) \
            <(echo -e HG01952#2#JAHAMD010000259.1"\t"chm13#chrX) \
            <(echo -e HG01978#2#JAGYVR010000154.1"\t"chm13#chrX) \
            <(echo -e HG02055#1#JAHEPK010000168.1"\t"chm13#chrY) \
            <(echo -e HG02055#1#JAHEPK010000257.1"\t"chm13#chrY) \
            <(echo -e HG02055#1#JAHEPK010000314.1"\t"chm13#chrY) \
            <(echo -e HG02145#1#JAHKSG010000262.1"\t"chm13#chrY) \
            <(echo -e HG02145#1#JAHKSG010000342.1"\t"chm13#chrY) \
            <(echo -e HG02145#1#JAHKSG010000399.1"\t"chm13#chrY) \
            <(echo -e HG02145#1#JAHKSG010000547.1"\t"chm13#chrY) \
            <(echo -e HG02145#1#JAHKSG010000739.1"\t"chm13#chrY) \
            <(echo -e HG02145#1#JAHKSG010000785.1"\t"chm13#chrY) \
            <(echo -e HG02257#1#JAGYVI010000079.1"\t"chm13#chrX) \
            <(echo -e HG02486#1#JAGYVM010000128.1"\t"chm13#chrY) \
            <(echo -e HG02486#1#JAGYVM010000157.1"\t"chm13#chrY) \
            <(echo -e HG02486#1#JAGYVM010000169.1"\t"chm13#chrY) \
            <(echo -e HG02486#1#JAGYVM010000291.1"\t"chm13#chrY) \
            <(echo -e HG02572#1#JAHAOW010000392.1"\t"chm13#chrY) \
            <(echo -e HG02572#1#JAHAOW010000399.1"\t"chm13#chrY) \
            <(echo -e HG02572#1#JAHAOW010000420.1"\t"chm13#chrY) \
            <(echo -e HG02572#1#JAHAOW010000474.1"\t"chm13#chrY) \
            <(echo -e HG02572#1#JAHAOW010000567.1"\t"chm13#chrY) \
            <(echo -e HG02717#1#JAHAOS010000264.1"\t"chm13#chrY) \
            <(echo -e HG02717#1#JAHAOS010000281.1"\t"chm13#chrY) \
            <(echo -e HG02717#1#JAHAOS010000282.1"\t"chm13#chrY) \
            <(echo -e HG02723#1#JAHEOU010000219.1"\t"chm13#chrX) \
            <(echo -e HG02886#2#JAHAOT010000399.1"\t"chm13#chrX) \
            <(echo -e HG03098#1#JAHEPM010000164.1"\t"chm13#chrY) \
            <(echo -e HG03098#1#JAHEPM010000254.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000041.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000272.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000452.1"\t"chm13#chrY) \
            <(echo -e HG03492#1#JAHEPI010000487.1"\t"chm13#chrY) \
            <(echo -e HG03579#1#JAGYVU010000262.1"\t"chm13#chrY) \
            <(echo -e HG03579#1#JAGYVU010000323.1"\t"chm13#chrY) \
            <(echo -e HG03579#1#JAGYVU010000368.1"\t"chm13#chrY) \
            <(echo -e HG03579#1#JAGYVU010000373.1"\t"chm13#chrY) \
            <(echo -e HG03579#1#JAGYVU010000417.1"\t"chm13#chrY) \
            <(echo -e HG01123#1#JAGYYZ010000259.1"\t"chm13#chrX) \
            <(echo -e HG002-bakeoff#PAT#scaffold_235"\t"chm13#chrY) \
            <(echo -e HG002-bakeoff#PAT#SY_unloc_11"\t"chm13#chrY) \
            <(echo -e HG002-bakeoff#PAT#SY_unloc_2"\t"chm13#chrY) \
            <(echo -e chm13#chrX"\t"chm13#chrX) \
            <(echo -e chm13#chrY"\t"chm13#chrY) \
            <(echo -e grch38#chrX"\t"chm13#chrX) \
            <(echo -e grch38#chrY"\t"chm13#chrY) > /lizardfs/guarracino/chromosome_communities/untangle_sex/partitioning.contig2chr.tsv


for e in 50000; do
  for m in 1000 ; do
    path_ref_fixed_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle_sex/$prefix.untangle.chm13#SEX.e$e.m$m.j0.n100.fixed.bed.gz
    if [[ ! -s ${path_ref_fixed_bed_gz} ]]; then
      echo "-e $e -m $m"
      
      path_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle_sex/$prefix.untangle.chm13#SEX.e$e.m$m.j0.n100.bed.gz
      python3 /lizardfs/guarracino/chromosome_communities/scripts/fix_best_hit.py \
        $path_bed_gz \
        /lizardfs/guarracino/chromosome_communities/untangle_sex/partitioning.contig2chr.tsv | tr ' ' '\t' | pigz -c -9 > $path_ref_fixed_bed_gz
    fi;
  done
done

# -p 98: Fixed 13900 hits covering 21791399 bps on the queries.
# -p 90: Fixed 12380 hits covering 20549819 bps on the queries.
```


Untangle with respect to a single sex chromosome, using always the same cut points:

```shell
for e in 50000; do
  for m in 1000; do
    echo "-e $e -m $m"
      
    cat $path_targets_txt | while read ref; do
      echo $ref
      
      path_ref_bed_gz=/lizardfs/guarracino/chromosome_communities/untangle_sex/$prefix.untangle.$ref.e$e.m$m.j0.n100.bed.gz
      if [[ ! -s ${path_ref_bed_gz} ]]; then
        path_cut_points_txt=/lizardfs/guarracino/chromosome_communities/untangle_sex/$prefix.untangle.chm13#SEX.e$e.m$m.j0.n100.cut_points.txt
        
        sbatch -p workers -c 48 --job-name sexuntangle --wrap "\time -v $RUN_ODGI untangle -t 48 -P -i $path_input_og -r $ref -e $e -m $m --cut-points-input $path_cut_points_txt -j 0 -n 100 | pigz -c > $path_ref_bed_gz"
      fi
    done
  done
done
```


Grounding (applying filters) and annotation:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded

for e in 50000; do
  for m in 1000; do
    cat $path_targets_txt | while read ref; do            
      path_grounded_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.$ref.e$e.m$m.grounded.tsv.gz
            
      if [[ ! -s ${path_grounded_tsv_gz} ]]; then
        echo "-e $e -m $m $ref grounding"

        # Grounding
        ( echo query query.begin query.end target target.begin target.end jaccard strand self.coverage nth.best ref ref.begin ref.end ref.jaccard ref.nth.best | tr ' ' '\t'
          join \
            <(zcat /lizardfs/guarracino/chromosome_communities/untangle_sex/$prefix.untangle.chm13#SEX.e$e.m$m.j0.n100.bed.gz | awk -v j=0 -v n=5 '{ if($7 >= j && $10 <= n) {print $1"_"$2, $0} }' | tr ' ' '\t' | sort -k 1,1) \
            <(zcat /lizardfs/guarracino/chromosome_communities/untangle_sex/$prefix.untangle.$ref.e$e.m$m.j0.n100.bed.gz | awk -v j=0 -v n=10 '{ if($7 >= j && $10 <= n) {print $1"_"$2, $0} }' | tr ' ' '\t' | sort -k 1,1) | \
          tr ' ' '\t' | grep -v '^#' | cut -f 2- | cut -f -10,14-17,20 | sort -k 1,3 -k 7,7nr -k 10,10n -k 14,14nr -k 15,15n ) | tr ' ' '\t' | pigz -c > x.tsv.gz


        echo "-e $e -m $m $ref filtering&annotation"
        ref_chr=$(echo $ref | sed 's/chm13#//')
        
        ###############################################################################
        # Contigs with respect to the partitioning
        zcat x.tsv.gz | sed '1d' | cut -f 1 | grep -v chr | grep -f <(grep $ref /lizardfs/guarracino/chromosome_communities/untangle_sex/partitioning.contig2chr.tsv | cut -f 1) | sort | uniq > $ref.tmp.txt

        # WRONG (it takes too many contigs): contigs overlapping (or close at least 100kbps to) a PARs/XTRs region
        #bedtools intersect \
        #  -a <(zcat x.tsv.gz | awk -v OFS="\t" '{print $11,$12,$13,$1, "", "+"}' | sed '1d' | bedtools sort) \
        #  -b <(grep $ref_chr /lizardfs/guarracino/chromosome_communities/data/chm13_hg002.PARs.bed |\
        #    bedtools sort |\
        #    bedtools slop -b 100000 -g /lizardfs/guarracino/chromosome_communities/data/chm13_hg002.PARs.sizes |\
        #    awk -v OFS='\t' -v ref=$ref '{print(ref,$2,$3)}') | \
        #  #awk '$3-$2+1>=100000' | \
        #  cut -f 4 | \
        #  #Remove references to avoid grepping everything later (with zgrep -f)
        #  grep -v chr |\
        #  sort | uniq > $ref.tmp.txt
        
        # WRONG (it takes too many contigs): contigs touching regions that are no PARs/XTRs
        #bedtools intersect \
        #  -a <(zcat x.tsv.gz | awk -v OFS="\t" '{print $11,$12,$13,$1, "", "+"}' | sed '1d' | grep "^$ref" | bedtools sort) \
        #  -b <(bedtools complement \
        #        -i <(grep $ref_chr data/chm13_hg002.PARs.bed | sed "s/$ref_chr/chm13/" | sed -e "s/PAR1/$ref_chr/" -e "s/PAR2/$ref_chr/" -e "s/XTR[12]/$ref_chr/" -e "s/XTR/$ref_chr/" | bedtools sort) -g <(cut -f 1,2 /lizardfs/guarracino/chromosome_communities/assemblies/chrSEX+refs.fa.gz.fai | grep $ref)            
        #  ) | \
        #  cut -f 4 | \
        #  #Remove references to avoid grepping everything later (with zgrep -f)
        #  grep -v chr |\
        #  sort | uniq > $ref.tmp.txt

        # WRONG: it mixes chrX and chrY contigs
        #zcat x.tsv.gz | sed '1d' | cut -f 1 | grep -v chr | sort | uniq > $ref.tmp.txt      
        ###############################################################################

        # Add grounded.target column, re-add the references, and add annotation
        cat \
          <(zcat x.tsv.gz | head -n 1 | awk -v OFS='\t' '{print $0, "grounded.target"}') \
          <(zgrep x.tsv.gz -f $ref.tmp.txt | awk -v OFS='\t' -v ref=$ref '{print $0, ref}') \
          <(zcat x.tsv.gz | awk -v OFS='\t' -v ref=$ref '$1 ~ /chr/ { print $0, ref}' ) \
          <(grep $ref_chr /lizardfs/guarracino/chromosome_communities/data/chm13_hg002.PARs.bed | \
            awk -v OFS='\t' -v ref=$ref '{print "chm13#"$1,".",".",ref,".",".","1","+","1","1",ref,$2,$3,"1","1",ref}') |\
          pigz -c -9 > $path_grounded_tsv_gz
        rm x.tsv.gz
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
               path_grounded_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.$ref.e$e.m$m.grounded.tsv.gz
      path_grounded_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.$ref.e$e.m$m.grounded.reliable.tsv.gz

      if [[ ! -s ${path_grounded_reliable_tsv_gz} ]]; then
        echo "-e $e -m $m $ref"

        # Skip verkko's and bakeoff's contigs because we don't have the unreliable regions for them
        zcat $path_grounded_tsv_gz | sed '1d' | grep 'HG002#MAT\|HG002#PAT\|HG01978#MAT\|HG01978#PAT\|bakeoff' -v | cut -f 1 | sort | uniq | while read CONTIG; do
          SAMPLE=$( echo $CONTIG | cut -f 1 -d '#')

          path_unreliable_bed=/lizardfs/guarracino/HPRC/annotations/unreliable/$SAMPLE.hifi.flagger_final.simplified.unreliable_only.bed
          if [[ -s $path_unreliable_bed ]]; then
            #echo $CONTIG "--->" $SAMPLE
            
            zgrep "^$CONTIG" $path_grounded_tsv_gz | sort -k 2n | awk -v OFS='\t' '{print($1,$2,$3,$4"_"$5"_"$6"_"$7"_"$8"_"$9"_"$10"_"$11"_"$12"_"$13"_"$14"_"$15"_"$16)}' > x.bed
            grep $CONTIG $path_unreliable_bed > y.bed
            # -A: remove entire feature if any overlap
            bedtools subtract -a x.bed -b y.bed -A |\
              awk -v OFS='\t' '{split($4, a, "_"); print($1,$2,$3,a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],a[11],a[12],a[13])}' >> x.tsv
            rm x.bed y.bed
          else
            zgrep "^$CONTIG" $path_grounded_tsv_gz | sort -k 2n >> x.tsv
          fi
        done
        
        # Re-take verkko's and bakeoff's untangled regions
        zcat $path_grounded_tsv_gz | sed '1d' | grep 'HG002#MAT\|HG002#PAT\|HG01978#MAT\|HG01978#PAT\|bakeoff' >> x.tsv
        
        cat \
          <( zcat $path_grounded_tsv_gz | head -n 1 ) \
          <( sort -k 1,3 -k 7,7nr -k 10,10n -k 14,14nr -k 15,15n x.tsv) | pigz -c -9 > $path_grounded_reliable_tsv_gz
        rm x.tsv
      fi;
    done
  done
done
```


Plot (`[start-500kbps,end+500kbps]` centered in the PARs/XTRs regions):

```shell
# PAR1/2/3
# https://link.springer.com/article/10.1007/s10142-013-0323-6/figures/1

n=1

for e in 50000; do
  for m in 1000; do
    for refn in 1; do
      (echo X; echo Y) | while read i; do     
        ref=chm13#chr$i
        
        echo "-e $e -m $m -refn $refn $ref"
    
        path_grounded_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.$ref.e$e.m$m.grounded.reliable.tsv.gz
        PREFIX=$(basename $path_grounded_reliable_tsv_gz .tsv.gz);
        
        zgrep '^HG01978#MAT\|^HG01978#PAT\|bakeoff' $path_grounded_reliable_tsv_gz -v | sed '1d' | cut -f 1 | sort | uniq \
          > /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.query_to_consider.txt
        
        if [[ $i == "X" ]]; then
            # chrX#PAR1:0-2394410
            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_without_annotation.R \
              $path_grounded_reliable_tsv_gz \
              0 2894410 \
              90 0.4 \
              0 \
              $n $refn \
              $i \
              0.9 \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.query_to_consider.txt \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.n$n.nref${refn}.PAR1.pdf
              
            # chrX#PAR2:153925834-154259566
            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_without_annotation.R \
              $path_grounded_reliable_tsv_gz \
              153425834 154259566 \
              90 0.4 \
              0 \
              $n $refn \
              $i \
              0.9 \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.query_to_consider.txt \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.n$n.nref${refn}.PAR2.pdf
              
            # chrX#PAR2:87642550-91570785
            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_without_annotation.R \
              $path_grounded_reliable_tsv_gz \
              87142550 92070785 \
              90 0.4 \
              0 \
              $n $refn \
              $i \
              0.9 \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.query_to_consider.txt \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.n$n.nref${refn}.XTR.pdf
        
            # Full chromosome X
            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_without_annotation.R \
              $path_grounded_reliable_tsv_gz \
              0 154259566 \
              200 0.4 \
              0 \
              $n $refn \
              $i \
              0.9 \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.query_to_consider.txt \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.n$n.nref${refn}.pdf
        else
            # chrY#PAR1:0-2458320
            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_without_annotation.R \
              $path_grounded_reliable_tsv_gz \
              0 2958320 \
              90 0.4 \
              0 \
              $n $refn \
              $i \
              0.9 \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.query_to_consider.txt \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.n$n.nref${refn}.PAR1.pdf
              
            # chrY#PAR2:62122809-62460029
            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_without_annotation.R \
              $path_grounded_reliable_tsv_gz \
              61622809 62460029 \
              90 0.4 \
              0 \
              $n $refn \
              $i \
              0.9 \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.query_to_consider.txt \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.n$n.nref${refn}.PAR2.pdf
              
            # chrY#XTR1:2727072-5914561
            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_without_annotation.R \
              $path_grounded_reliable_tsv_gz \
              2227072 6414561 \
              90 0.4 \
              0 \
              $n $refn \
              $i \
              0.9 \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.query_to_consider.txt \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.n$n.nref${refn}.XTR1.pdf
              
            # chrY#XTR2:6200973-6400875
            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_without_annotation.R \
              $path_grounded_reliable_tsv_gz \
              5700973 6900875 \
              90 0.4 \
              0 \
              $n $refn \
              $i \
              0.9 \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.query_to_consider.txt \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.n$n.nref${refn}.XTR2.pdf
            
            # Full chromosome Y
            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_without_annotation.R \
              $path_grounded_reliable_tsv_gz \
              0 62460029 \
              80 0.4 \
              0 \
              $n $refn \
              $i \
              0.9 \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.query_to_consider.txt \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.n$n.nref${refn}.pdf
        fi
      done
      
      # Merge chromosomes's PDF files
      #rm /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.chm13#chrSEX.e$e.m$m.grounded.reliable.n$n.nref${refn}.merged.pdf
      #/gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite \
      #  /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.chm13#chr*.e$e.m$m.grounded.reliable.n$n.nref${refn}.*pdf \
      #  /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.chm13#chrSEX.e$e.m$m.grounded.reliable.n$n.nref${refn}.merged.pdf
      #rm /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.chm13#chr*.e$e.m$m.grounded.reliable.n$n.nref${refn}.pdf
      #rm /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.chm13#chr*.e$e.m$m.grounded.reliable.n$n.nref${refn}.*R*.pdf
    done
  done
done


# Full chromosome X
for e in 50000; do
  for m in 1000; do
    for refn in 1; do
      (echo X; echo Y) | while read i; do     
        ref=chm13#chr$i
        
        echo "-e $e -m $m -refn $refn $ref"
    
        path_grounded_reliable_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.$ref.e$e.m$m.grounded.reliable.tsv.gz
        PREFIX=$(basename $path_grounded_reliable_tsv_gz .tsv.gz);
        
        zgrep '^HG01978#MAT\|^HG01978#PAT\|bakeoff' $path_grounded_reliable_tsv_gz -v | sed '1d' | cut -f 1 | sort | uniq \
          > /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.query_to_consider.txt
        
        if [[ $i == "X" ]]; then        
            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_without_annotation.R \
              $path_grounded_reliable_tsv_gz \
              0 154259566 \
              200 0.4 \
              0 \
              2 $refn \
              $i \
              0.7 \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.query_to_consider.txt \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.n2.nref${refn}.eid07.pdf
        else           
            Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_untangle_without_annotation.R \
              $path_grounded_reliable_tsv_gz \
              0 62460029 \
              80 0.4 \
              0 \
              2 $refn \
              $i \
              0.7 \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.query_to_consider.txt \
              /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$PREFIX.n2.nref${refn}.eid07.pdf
        fi
      done
    done
  done
done
```

CONTINUE (IF NEEDED) - Plot 2 hits: XXX

CONTINUE (IF NEEDED) - Compute support: XXX

```shell
# Merge files for all sex chromosomes (used for computing the support and the histogram length)
for e in 50000; do
  for m in 1000 ; do
    path_grounded_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.ALL.e$e.m$m.grounded.reliable.tsv.gz
    if [[ ! -s ${path_grounded_reliable_ALL_tsv_gz} ]]; then
      cat \
        <(zcat /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.chm13#chr*.e$e.m$m.grounded.reliable.tsv.gz | grep 'self.coverage' -m 1) \
        <(zcat /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.chm13#chr*.e$e.m$m.grounded.reliable.tsv.gz | grep 'self.coverage' -v ) |\
        pigz -c -9 > $path_grounded_reliable_ALL_tsv_gz
    fi;
  done
done

#Take sex chromosome lengths
grep '^chm13' /lizardfs/guarracino/chromosome_communities/assemblies/chm13v2+grch38masked.fa.gz.fai | grep 'chrX\|chrY' | cut -f 1,2 \
  > /lizardfs/guarracino/chromosome_communities/chm13#SEX.len.tsv

# Support
# guix install r-ggridges
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx todo
```


CONTINUE (IF NEEDED) - Statistics on removed regions: XXX


CONTINUE (IF NEEDED) - Statistics on untangled segment lengths: XXX

TO FINISH, IF NEEDED
Estimate regions that can recombine using multi-hit untangled regions:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/recombinant_regions/

# Identify regions with multiple good enough hits
for e in 50000; do
  for m in 1000; do
    path_grounded_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.ALL.e$e.m$m.grounded.reliable.tsv.gz
    PREFIX=$(basename $path_grounded_reliable_ALL_tsv_gz .tsv.gz)
   
    zgrep '^chm13\|^grch38\|^HG002#1\|HG002#2\|^HG01978#MAT\|^HG01978#PAT\|bakeoff' $path_grounded_reliable_ALL_tsv_gz -v | sed '1d' | cut -f 1 | sort | uniq \
      > /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/recombinant_regions/$PREFIX.query_to_consider.txt

    for sc in 0 1.1; do
      for eid in 0.900 0.950 0.975 0.995 1.000; do
        eid_str=$(echo $eid | sed 's/\.//g')
        sc_str=$(echo $sc | sed 's/\.//g')
        echo $e $m $sc $eid
          
        path_recombinant_regions_bed=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/recombinant_regions/$PREFIX.recombinant_regions.sc${sc_str}.eid${eid_str}.bed
        if [[ ! -s $path_recombinant_regions_bed ]]; then
          python3 /lizardfs/guarracino/chromosome_communities/scripts/recombination_proxy_ranges.py \
            $path_grounded_reliable_ALL_tsv_gz \
            $eid $sc \
            /lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/recombinant_regions/$PREFIX.query_to_consider.txt \
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
awk -v OFS='\t' '{print($1,"0",$2)}' /lizardfs/guarracino/chromosome_communities/chm13#SEX.len.tsv > chm13.bed

path_targets_txt=/lizardfs/guarracino/chromosome_communities/untangle_sex/chm13.SEX.target_paths.txt
path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrSEX+refs.s50k.l250k.p98.n102/chrSEX+refs.fa.gz.2ed2c67.04f1c29.22fc5c8.smooth.final.og
prefix=$(basename $path_input_og .og)

for e in 50000; do
  for m in 1000; do
    path_grounded_reliable_ALL_tsv_gz=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/$prefix.untangle.ALL.e$e.m$m.grounded.reliable.tsv.gz
    PREFIX=$(basename $path_grounded_reliable_ALL_tsv_gz .tsv.gz)

    path_recombinant_regions_table_tsv=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/recombinant_regions/$PREFIX.recombinant_regions.table.tsv
    path_recombinant_regions_table_sizes_tsv=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/recombinant_regions/$PREFIX.recombinant_regions.table.sizes.tsv
    path_recombinant_regions_table_with_counts_tsv=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/recombinant_regions/$PREFIX.recombinant_regions.table.counts.tsv
    
    rm rg.txt
    rm $path_recombinant_regions_table_tsv $path_recombinant_regions_table_sizes_tsv $path_recombinant_regions_table_with_counts_tsv
    for sc in 0 1.1; do
      for eid in 0.900 0.950 0.975 0.995 1.000; do
        eid_str=$(echo $eid | sed 's/\.//g')
        sc_str=$(echo $sc | sed 's/\.//g')
        echo $e $m $sc $eid
          
        path_recombinant_regions_bed=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/recombinant_regions/$PREFIX.recombinant_regions.sc${sc_str}.eid${eid_str}.bed

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
# Use the average counts as self coverage for the putative recombinant regions
cd /lizardfs/guarracino/chromosome_communities/
f=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/recombinant_regions/chrSEX+refs.fa.gz.2ed2c67.04f1c29.22fc5c8.smooth.final.untangle.ALL.e50000.m1000.grounded.reliable.recombinant_regions.table.counts.tsv
max_count=$(cut -f 6 $f | sort -n | tail -n 1)
echo -e "query\tquery.begin\tquery.end\ttarget\ttarget.begin\ttarget.end\tjaccard\tstrand\tself.coverage\tnth.best\tref\tref.begin\tref.end\tref.jaccard\tref.nth.best\tgrounded.target" > chrXY+recombinant.tsv
for chr in chrX chrY; do
  cat \
    <( grep $chr /lizardfs/guarracino/chromosome_communities/data/chm13_hg002.PARs.bed | sed 's/^/chm13#/g' | \
        awk -v OFS='\t' -v ref=chm13#$chr -v max=$max_count '{print $1,".",".",ref,".",".","1","+",max,"1",ref,$2,$3,"1","1",ref}' ) \
    <( grep $chr $f | \
        awk '$1 == 0 && $2 == 0.9' | cut -f 3,4,5,6 | awk '$4 > 0' | bedtools sort | bedtools merge -o mean -c 4 -d 10000 | \
        awk '$3 - $2 > 30000' | \
        awk -v OFS='\t' -v ref=chm13#$chr '{print $1,".",".",ref,".",".","1","+",$4,"1",ref,$2,$3,"1","1",ref}' ) | sort -r >> chrXY+recombinant.tsv
done

Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_recombinant_regions_chrXY.R \
  chrXY+recombinant.tsv \
  90 10 \
  "X" \
  /lizardfs/guarracino/chromosome_communities/chrX.recombinant_regions.eid0900.pdf # Supplementary figure
Rscript /lizardfs/guarracino/chromosome_communities/scripts/plot_recombinant_regions_chrXY.R \
  chrXY+recombinant.tsv \
  90 12 \
  "Y" \
  /lizardfs/guarracino/chromosome_communities/chrY.recombinant_regions.eid0900.pdf # Supplementary figure
  
# Merge chromosomes's PDF files
/gnu/store/d0njxcgymxvf8s7di32m9q4v9vibd11z-poppler-0.86.1/bin/pdfunite \
  /lizardfs/guarracino/chromosome_communities/chrX.recombinant_regions.eid0900.pdf /lizardfs/guarracino/chromosome_communities/chrY.recombinant_regions.eid0900.pdf \
  /lizardfs/guarracino/chromosome_communities/chrXY.recombinant_regions.eid0900.pdf
rm /lizardfs/guarracino/chromosome_communities/chrX.recombinant_regions.eid0900.pdf /lizardfs/guarracino/chromosome_communities/chrY.recombinant_regions.eid0900.pdf



f=/lizardfs/guarracino/chromosome_communities/untangle_sex/grounded/recombinant_regions/chrSEX+refs.fa.gz.2ed2c67.04f1c29.22fc5c8.smooth.final.untangle.ALL.e50000.m1000.grounded.reliable.recombinant_regions.table.counts.tsv
awk '$1 == 0 && $2 == 0.9 && $6 > 0' $f | cut -f 3,4,5,6 | awk '$4 > 0' | bedtools sort | bedtools merge | awk -v OFS='\t' '{SUM+=$3-$2}END{print(SUM)}'
awk '$1 == 0 && $2 == 0.9 && $6 > 0' $f | cut -f 3,4,5,6 | awk '$4 > 0' | bedtools sort | bedtools merge | grep chrX | awk -v OFS='\t' '{SUM+=$3-$2}END{print(SUM)}'
awk '$1 == 0 && $2 == 0.9 && $6 > 0' $f | cut -f 3,4,5,6 | awk '$4 > 0' | bedtools sort | bedtools merge | grep chrY | awk -v OFS='\t' '{SUM+=$3-$2}END{print(SUM)}'

awk '$1 == 0 && $2 == 0.9 && $6 > 19' $f | cut -f 3,4,5,6 | awk '$4 > 0' | bedtools sort | bedtools merge | awk -v OFS='\t' '{SUM+=$3-$2}END{print(SUM)}'
awk '$1 == 0 && $2 == 0.9 && $6 > 19' $f | cut -f 3,4,5,6 | awk '$4 > 0' | bedtools sort | bedtools merge | grep chrX | awk -v OFS='\t' '{SUM+=$3-$2}END{print(SUM)}'
awk '$1 == 0 && $2 == 0.9 && $6 > 19' $f | cut -f 3,4,5,6 | awk '$4 > 0' | bedtools sort | bedtools merge | grep chrY | awk -v OFS='\t' '{SUM+=$3-$2}END{print(SUM)}'
```


[//]: # (## Variant calling)
[//]: # ()
[//]: # (Call variants in a haploid setting:)
[//]: # ()
[//]: # (```shell)
[//]: # (path_input_gfa=/lizardfs/guarracino/chromosome_communities/graphs/chrS.pan+HG002chrXY.s100k.l300k.p98.n93/chrS.pan+HG002chrXY.fa.gz.73e7992.4030258.b8e2fe5.smooth.fix.gfa)
[//]: # (path_chm13_vcf_gz=/lizardfs/guarracino/chromosome_communities/graphs/chrS.pan+HG002chrXY.s100k.l300k.p98.n93/chrS.pan+HG002chrXY.fa.gz.73e7992.4030258.b8e2fe5.smooth.fix.chm13.vcf.gz)
[//]: # (#path_grch38_vcf_gz=/lizardfs/guarracino/chromosome_communities/graphs/chrS.pan+HG002chrXY.s100k.l300k.p98.n93/chrS.pan+HG002chrXY.fa.gz.73e7992.4030258.b8e2fe5.smooth.fix.grch38.vcf.gz)
[//]: # (path_hg002_vcf_gz=/lizardfs/guarracino/chromosome_communities/graphs/chrS.pan+HG002chrXY.s100k.l300k.p98.n93/chrS.pan+HG002chrXY.fa.gz.73e7992.4030258.b8e2fe5.smooth.fix.hg002.vcf.gz)
[//]: # (sbatch -p workers -c 48 --job-name vgchm13 --wrap '\time -v vg deconstruct -P chm13 -H '?' -e -a -t 48 '$path_input_gfa' | bgzip -@ 48 -c > '$path_chm13_vcf_gz' && tabix '$path_chm13_vcf_gz)
[//]: # (#sbatch -p workers -c 48 --job-name vggrch38 --wrap '\time -v vg deconstruct -P grch38 -H '?' -e -a -t 48 '$path_input_gfa' | bgzip -@ 48 -c > '$path_grch38_vcf_gz' && tabix '$path_grch38_vcf_gz)
[//]: # (sbatch -p workers -c 48 --job-name vghg002 --wrap '\time -v vg deconstruct -P HG002 -H '?' -e -a -t 48 '$path_input_gfa' | bgzip -@ 48 -c > '$path_hg002_vcf_gz' && tabix '$path_hg002_vcf_gz)
[//]: # ()
[//]: # ()
[//]: # ()
[//]: # ()
[//]: # ()
[//]: # (# In the VCF there are variants with all samples having missing genotype!)
[//]: # (#num_samples=`bcftools query -l $PATH_VCF_GZ | wc -l`)
[//]: # (#num_miss_gen=$&#40;echo $num_samples - 1 | bc&#41;)
[//]: # (#--max-missing-count $num_miss_gen \)
[//]: # (```)
