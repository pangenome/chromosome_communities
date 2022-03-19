### msmc

```shell

cd ~/tools
git clone --recursive https://github.com/stschiff/msmc-tools.git

cd ~/tools
wget -c https://github.com/stschiff/msmc/releases/download/v1.1.0/msmc_1.1.0_linux64bit
chmod +x msmc_1.1.0_linux64bit
```

Tutorial: https://github.com/stschiff/msmc-tools/blob/master/msmc-tutorial/guide.md#estimating-population-separation-history


```shell
#grep -e $'chr1\t\|chr2\t\|chr3\|chr4\|chr5\|cht6\|cht7\|chr8\|chr9\|chr13\|chr14\|chr15\|chr21\|chr22'
THREADS=48

PATH_SNV_VCF_GZ=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.snv.vcf.gz
VCF_NAME=$(basename $PATH_SNV_VCF_GZ .vcf.gz)

mkdir -p /lizardfs/guarracino/chromosome_communities/msmc
cd /lizardfs/guarracino/chromosome_communities/msmc

# Split p/q-arm variants 
sed 's/chr/chm13#chr/g' /lizardfs/guarracino/chromosome_communities/data/chm13.centromeres.approximate.bed | \
  bedtools sort | \
  bedtools complement \
    -i - \
    -g <(cut -f 1,2 /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa.fai | sort) | \
  grep 'chr13\|chr14\|chr15\|chr21\|chr22' > tmp.bed
  
# NOTE: in vctools, the BED file is expected to have a header line.
# Take odd rows
(echo -e "#chrom\tstart\tend"; sed -n 1~2p tmp.bed) > /lizardfs/guarracino/chromosome_communities/msmc/p_arms.bed
# Take even rows
(echo -e "#chrom\tstart\tend"; sed -n 2~2p tmp.bed) > /lizardfs/guarracino/chromosome_communities/msmc/q_arms.bed

rm tmp.bed


vcftools \
  --gzvcf $PATH_SNV_VCF_GZ \
  --bed /lizardfs/guarracino/chromosome_communities/msmc/p_arms.bed \
  --recode --out p --keep-INFO-all --stdout |\
  #bcftools norm -f $PATH_REF_FASTA -c e -m - |\
  bgzip -@ $THREADS > $VCF_NAME.p_arms.vcf.gz
vcftools \
  --gzvcf $PATH_SNV_VCF_GZ \
  --bed /lizardfs/guarracino/chromosome_communities/msmc/q_arms.bed \
  --recode --out q --keep-INFO-all --stdout |\
  #bcftools norm -f $PATH_REF_FASTA -c e -m - |\
  bgzip -@ $THREADS > $VCF_NAME.q_arms.vcf.gz
rm *log

# Split by sample and fill missing genotypes (with 0, i.e. reference allele, else msmc breaks)
for SAMPLE in `bcftools query -l $PATH_SNV_VCF_GZ`; do
    for ARMS in p_arms q_arms; do
      echo $SAMPLE $ARMS
      bcftools view --samples $SAMPLE $VCF_NAME.$ARMS.vcf.gz |\
        bcftools +setGT - -- -t . --new-gt 0p |\
        bgzip -@ 48 -c > $VCF_NAME.$ARMS.$SAMPLE.vcf.gz
    done
done

# Take samples from each ACRO
rm samples.p_arms.txt samples.q_arms.txt

#nsamples=15   # Take 15 samples
#(seq 13 15; seq 21 22) | while read i; do
#  PATH_REF_FASTA_FAI=/lizardfs/guarracino/chromosome_communities/pq_contigs/chr$i.vs.chm13.50kbps.pq_contigs.union.hg002prox.fa.gz.fai
#
#  head $PATH_REF_FASTA_FAI -n $nsamples | cut -f 1 | while read SAMPLE; do
#    echo "chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.snv.p_arms.${SAMPLE}.vcf.gz" >> samples.p_arms.txt
#    echo "chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.snv.q_arms.${SAMPLE}.vcf.gz" >> samples.q_arms.txt
#  done
#done

(seq 13 15; seq 21 22) | while read i; do
  PATH_REF_FASTA_FAI=/lizardfs/guarracino/chromosome_communities/pq_contigs/chr$i.vs.chm13.50kbps.pq_contigs.union.hg002prox.fa.gz.fai

  cat $PATH_REF_FASTA_FAI $nsamples | cut -f 1 | while read SAMPLE; do
    echo "chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.snv.p_arms.${SAMPLE}.vcf.gz" >> samples.p_arms.txt
    echo "chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.snv.q_arms.${SAMPLE}.vcf.gz" >> samples.q_arms.txt
  done
done
nsamples=$(cat samples.p_arms.txt | wc -l)
nsamples_minus_1=$(echo "$nsamples - 1" | bc)
nsamples_mul_2_minus_1=$(echo "$nsamples * 2 - 1" | bc)


samples=$(cat samples.p_arms.txt samples.p_arms.txt | tr '\n' ' ' | sed 's/.$//'`)

# Generate files for MSMC
# Using a hacked version, changing a few rows in generate_multihetsep.py to get phased fake-diploid data
#print ("HACKED: Non-diploid SNP found and considered as phased data: %s" % geno, file=sys.stderr)
#phased = True
#geno = "%s|%s" % (geno[0], geno[0])
python3 ~/tools/msmc-tools/generate_multihetsep.py $samples > $nsamples.pq_arms.phased.hetsep

# msmc wants a single chromosome
# sed to avoid '#' in the reference name that breaks msmc
grep '^chm13#chr13' $nsamples.pq_arms.phased.hetsep | sed 's/^chm13#//g' > $nsamples.pq_arms.phased.chr13.hetsep

#!/bin/sh

# Estimating the effective population size
haplos1=`seq 0 $nsamples_minus_1 | tr '\n' ',' | sed 's/.$//'`
#sbatch -p allnodes -c 48 --job-name msmc2_p1 --wrap 'hostname; \time -v /home/guarracino/tools/msmc2_linux64bit -p 1*2+15*1+1*2 -o '$nsamples'.pq_arms.phased.chr13.hetsep.p_arms -I '$haplos1' /lizardfs/guarracino/chromosome_communities/msmc/'$nsamples'.pq_arms.phased.chr13.hetsep'
\time -v /home/guarracino/tools/msmc2_linux64bit -p 1*2+15*1+1*2 -o $nsamples.pq_arms.phased.chr13.hetsep.p_arms -I $haplos1 /lizardfs/guarracino/chromosome_communities/msmc/$nsamples.pq_arms.phased.chr13.hetsep

haplos2=`seq $nsamples $nsamples_mul_2_minus_1 | tr '\n' ',' | sed 's/.$//'`
#sbatch -p allnodes -c 48 --job-name msmc2_p2 --wrap 'hostname; \time -v /home/guarracino/tools/msmc2_linux64bit -p 1*2+15*1+1*2 -o '$nsamples'.pq_arms.phased.chr13.hetsep.q_arms -I '$haplos2' /lizardfs/guarracino/chromosome_communities/msmc/'$nsamples'.pq_arms.phased.chr13.hetsep'
\time -v /home/guarracino/tools/msmc2_linux64bit -p 1*2+15*1+1*2 -o $nsamples.pq_arms.phased.chr13.hetsep.q_arms -I $haplos2 /lizardfs/guarracino/chromosome_communities/msmc/$nsamples.pq_arms.phased.chr13.hetsep


# Estimating population separation history
combinations=$((for x in `seq 0 $nsamples_minus_1`; do
  for y in `seq $nsamples $nsamples_mul_2_minus_1`; do
    echo "$x-$y"
  done
done) | tr '\n' ',' | sed 's/.$//')

# --skipAmbiguous: skip sites with ambiguous phasing. Recommended for cross population analysis
#sbatch -p allnodes -c 48 --job-name msmc2_sep12 --wrap 'hostname; \time -v ~/tools/msmc2_linux64bit -p 1*2+15*1+1*2 --skipAmbiguous -o '$nsamples'.pq_arms.phased.chr13.hetsep.pq_arms -I '$combinations' /lizardfs/guarracino/chromosome_communities/msmc/'$nsamples'.pq_arms.phased.chr13.hetsep'
\time -v ~/tools/msmc2_linux64bit -p 1*2+15*1+1*2 --skipAmbiguous -o $nsamples.pq_arms.phased.chr13.hetsep.pq_arms -I $combinations /lizardfs/guarracino/chromosome_communities/msmc/$nsamples.pq_arms.phased.chr13.hetsep


python3 ~/tools/msmc-tools/combineCrossCoal.py \
  $nsamples.pq_arms.phased.chr13.hetsep.pq_arms.final.txt \
  $nsamples.pq_arms.phased.chr13.hetsep.p_arms.final.txt \
  $nsamples.pq_arms.phased.chr13.hetsep.q_arms.final.txt \
  > $nsamples.pq_arms.phased.chr13.hetsep.combined.final.txt








# sed to avoid '#' in the reference name that breaks msmc
python3 ~/tools/msmc-tools/generate_multihetsep.py \
  chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.snv.HG002#MAT#chr13.prox.vcf.gz\
  chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.snv.HG00438#1#JAHBCB010000045.1.vcf.gz \
  chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.snv.HG00621#2#JAHBCC010000104.1.vcf.gz \
  chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.snv.HG00741#1#JAHALY010000001.1.vcf.gz \
  chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.snv.HG01978#1#JAGYVS010000056.1.vcf.gz \
  chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.snv.HG002#MAT#chr14.prox.vcf.gz\
  chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.snv.HG00673#2#JAHBBY010000065.1.vcf.gz \
  chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.snv.HG00733#1#JAHEPQ010000013.1.vcf.gz \
  chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.snv.HG00735#1#JAHBCH010000039.1.vcf.gz \
  chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.snv.HG01978#1#JAGYVS010000055.1.vcf.gz | \
  sed 's/^chm13#//g' > subset.13vs14.phased.hetsep

# msmc wants a single chromosome
grep 'chr13' subset.13vs14.phased.hetsep > subset.13vs14.phased.chr13.hetsep

# Estimating the effective population size
~/tools/msmc2_linux64bit -p 1*2+15*1+1*2 -o subset.13vs14.chr13.pop13 subset.13vs14.phased.chr13.hetsep -I 0,1,2,3,4 
~/tools/msmc2_linux64bit -p 1*2+15*1+1*2 -o subset.13vs14.chr13.pop14 subset.13vs14.phased.chr13.hetsep -I 5,6,7,8,9

# Estimating population separation history
~/tools/msmc2_linux64bit -p 1*2+15*1+1*2 -o subset.13vs14.chr13.separation subset.13vs14.phased.chr13.hetsep \
  -I 0-5,0-6,0-7,0-8,0-9,1-5,1-6,1-7,1-8,1-9,2-5,2-6,2-7,2-9,2-9,3-5,3-6,3-7,3-8,3-9,4-5,4-6,4-7,4-8,4-9


python3 ~/tools/msmc-tools/combineCrossCoal.py subset.13vs14.chr13.separation.final.txt \
  subset.13vs14.chr13.pop13.final.txt \
  subset.13vs14.chr13.pop14.final.txt > subset.13vs14.chr13.separation.combined.final.txt











python3 ~/tools/msmc-tools/generate_multihetsep.py \
  chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.chm13.*.vcf.gz \
  > sed 's/chm13#//g' chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.chm13.ALL.hetsep
# The first two columns denote chromosome and position of a segregating site within the samples.
# The fourth column contains the X alleles in the X(?) haplotypes we put in. 
# When there are multiple patterns separated by a comma, it means that phasing information is ambiguous, 
# so there are multiple possible phasings.


~/tools/msmc2_linux64bit -p 1*2+15*1+1*2 -o prefox.chr13 chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.chm13.ALL.hetsep -I 0,1




```