# History - PGGB + MSMC2/SMC++

## Tools and directory

```shell
cd ~/tools
git clone --recursive https://github.com/stschiff/msmc-tools.git
# Hack the tool, changing a few rows in generate_multihetsep.py to get phased fake-diploid data
#if len(geno) != 3 or geno[1] not in ['|/']:
#  print ("HACKED: Non-diploid SNP found and considered as phased data: %s" % geno, file=sys.stderr)
#  phased = True
#  geno = "%s|%s" % (geno[0], geno[0])
cd ..

wget -c https://github.com/stschiff/msmc/releases/download/v1.1.0/msmc_1.1.0_linux64bit
chmod +x msmc_1.1.0_linux64bit

wget -c https://github.com/popgenmethods/smcpp/releases/download/v1.15.2/smcpp-1.15.2-Linux-x86_64.sh
bash smcpp-1.15.2-Linux-x86_64.sh
/home/guarracino/smcpp/bin/smc++ 
```


## Notes
In vcftools, the BED file is expected to have a header line.

MSMC2 tutorial: https://github.com/stschiff/msmc-tools/blob/master/msmc-tutorial/guide.md#estimating-population-separation-history


## Preparation

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/msmc
cd /lizardfs/guarracino/chromosome_communities/msmc
```

Prepare p/q-arms BED files:

```shell
sed 's/chr/chm13#chr/g' /lizardfs/guarracino/chromosome_communities/data/chm13.centromeres.approximate.bed | \
  bedtools sort | \
  bedtools complement \
    -i - \
    -g <(cut -f 1,2 /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa.fai | grep 'chrM\|chrX' -v | sort) \
  > tmp.bed

# Take odd rows
(echo -e "#chrom\tstart\tend"; sed -n 1~2p tmp.bed | grep 'chr13\|chr14\|chr15\|chr21\|chr22') > p_arms.bed

# Take even rows
(echo -e "#chrom\tstart\tend"; sed -n 2~2p tmp.bed | grep 'chr13\|chr14\|chr15\|chr21\|chr22') > q_arms.bed
rm tmp.bed


# Regions on the right of the rDNA
zgrep rDNA /lizardfs/guarracino/chromosome_communities/data/chm13.CentromeresAndTelomeres.CenSatAnnotation.txt.gz | \
  sed 's/chr/chm13#chr/g' | \
  awk -v OFS='\t' '{print $1,"0",$3}' |\
  bedtools complement \
    -i - \
    -g <(grep 'chr13\|chr14\|chr15\|chr21\|chr22' /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa.fai | cut -f 1,2 | sort) \
  > rDNA.right.bed

bedtools intersect -a p_arms.bed -b rDNA.right.bed > p_arms.right.bed
```


## All possible pairs of acrocentric chromosomes

Prepare acro-VCF files for each sample, with respect to each acro-reference:

```shell
PATH_VCF_GZ=/lizardfs/guarracino/chromosome_communities/graphs/chrA.pan+HG002chrAprox.s100k.l300k.p98.n100/chrA.pan+HG002chrAprox.fa.gz.952ce3b.4030258.0e65f78.smooth.fix.sed.c1000.chm13.vcf.gz
VCF_NAME=$(basename $PATH_VCF_GZ .vcf.gz)

# Consider p-arm-right SNVs with at most 70% of missing genotypes
(seq 13 15; seq 21 22) | while read i; do
  REF=chr$i
  echo $REF

  mkdir -p $REF

  # Do not consider the GRCh38 reference (very bad in on the p-arm of the acrocentric chromosomes)
  vcftools \
    --gzvcf $PATH_VCF_GZ \
    --bed <(echo -e "#chrom\tstart\tend"; grep $REF q_arms.bed | sed 's/#/-/') \
    --remove-indv 'grch38-chr13' \
    --remove-indv 'grch38-chr14' \
    --remove-indv 'grch38-chr15' \
    --remove-indv 'grch38-chr21' \
    --remove-indv 'grch38-chr22' \
    --max-missing 1 \
    --remove-indels \
    --recode --keep-INFO-all \
    --out $REF.p_arm.right.maxMiss1 \
    --stdout |\
    #bcftools norm -f $PATH_REF_FASTA -c e -m - |\
    bgzip -@ 48 > $REF/$VCF_NAME.$REF.p_arm.right.vcf.gz 
done

#(seq 13 15; seq 21 22) | while read i; do
#  REF=chr$i
#  echo $REF
#
#  # Do not consider the GRCh38 reference (very bad in on the p-arm of the acrocentric chromosomes)
#  vcftools \
#    --gzvcf $PATH_VCF_GZ \
#    --bed <(echo -e "#chrom\tstart\tend"; grep $REF q_arms.bed | sed 's/#/-/') \
#    --remove-indv 'grch38-chr13' \
#    --remove-indv 'grch38-chr14' \
#    --remove-indv 'grch38-chr15' \
#    --remove-indv 'grch38-chr21' \
#    --remove-indv 'grch38-chr22' \
#    --remove-indels \
#    --recode --keep-INFO-all \
#    --stdout | bcftools stats -s '-' - > $REF.stats.tsv
#done

#--max-missing 1.0
#chr13.p_arm.right.maxMiss1.log:After filtering, kept 620046 out of a possible 3708776 Sites
#chr14.p_arm.right.maxMiss1.log:After filtering, kept 544422 out of a possible 3708776 Sites
#chr15.p_arm.right.maxMiss1.log:After filtering, kept 504982 out of a possible 3708776 Sites
#chr21.p_arm.right.maxMiss1.log:After filtering, kept 234702 out of a possible 3708776 Sites
#chr22.p_arm.right.maxMiss1.log:After filtering, kept 228922 out of a possible 3708776 Sites

#--max-missing 0.9
#chr13.p_arm.right.maxMiss09.log:After filtering, kept 645105 out of a possible 3708776 Sites
#chr14.p_arm.right.maxMiss09.log:After filtering, kept 584216 out of a possible 3708776 Sites
#chr15.p_arm.right.maxMiss09.log:After filtering, kept 538064 out of a possible 3708776 Sites
#chr21.p_arm.right.maxMiss09.log:After filtering, kept 246852 out of a possible 3708776 Sites
#chr22.p_arm.right.maxMiss09.log:After filtering, kept 256023 out of a possible 3708776 Sites

#--max-missing 0.0
#chr13.p_arm.right.maxMiss0.log:After filtering, kept 656414 out of a possible 3708776 Sites
#chr14.p_arm.right.maxMiss0.log:After filtering, kept 608339 out of a possible 3708776 Sites
#chr15.p_arm.right.maxMiss0.log:After filtering, kept 575662 out of a possible 3708776 Sites
#chr21.p_arm.right.maxMiss0.log:After filtering, kept 251152 out of a possible 3708776 Sites
#chr22.p_arm.right.maxMiss0.log:After filtering, kept 270264 out of a possible 3708776 Sites

```


```shell



# Split p/q-arm variants 

(seq 13 15; seq 21 22) | while read i; do
  REF=chr$i
  echo $REF

  mkdir -p $REF
    
  vcftools \
    --gzvcf $PATH_SNV_VCF_GZ \
    --bed <(echo -e "#chrom\tstart\tend"; grep $REF p_arms.bed) \
    --recode --keep-INFO-all --out $REF.p --stdout |\
    #bcftools norm -f $PATH_REF_FASTA -c e -m - |\
    bgzip -@ 48 > $REF/$VCF_NAME.$REF.p_arm.vcf.gz
  vcftools \
    --gzvcf $PATH_SNV_VCF_GZ \
    --bed <(echo -e "#chrom\tstart\tend"; grep $REF q_arms.bed) \
    --recode --keep-INFO-all --out $REF.q --stdout |\
    #bcftools norm -f $PATH_REF_FASTA -c e -m - |\
    bgzip -@ 48 > $REF/$VCF_NAME.$REF.q_arm.vcf.gz
done
#rm *log




# Split by sample and fill missing genotypes (with 0, i.e. reference allele, else msmc breaks)
(seq 13 15; seq 21 22) | while read i; do
  REF=chr$i
    
  for ARM in p_arm q_arm; do
    #Better to keep the same samples in p/q VCFs
    #SAMPLES_WITH_GENOTYPES=$(bcftools stats -s '-' $REF/$VCF_NAME.$REF.$ARM.vcf.gz | awk '$1=="PSC" && $12+$13>0 {print $3}')
    for SAMPLE in `bcftools query -l $PATH_SNV_VCF_GZ | tr '\n' ' '`; do
      echo $REF $SAMPLE $ARM
      bcftools view --samples $SAMPLE $REF/$VCF_NAME.$REF.$ARM.vcf.gz |\
        bcftools +setGT - -- -t . --new-gt 0p |\
        bgzip -@ 48 -c > $REF/$VCF_NAME.$REF.$ARM.$SAMPLE.vcf.gz
    done
  done
done
```









Split variants by arm:

```shell
(seq 1 22) | while read i; do
  REF=chr$i
  echo $REF 

  PATH_CHR_VCF_GZ=/lizardfs/erikg/HPRC/year1v2genbank/wgg.88/$REF.pan/$REF.pan.fa.a2fb268.4030258.6a1ecc2.smooth.chm13.vcf.gz
  VCF_NAME=$(basename $PATH_CHR_VCF_GZ .vcf.gz)

  for ARM in p_arm q_arm; do
    ARM_DIR=$REF/$ARM
    mkdir -p $ARM_DIR

    PATH_CHR_ARM_SNV_VCF_GZ=$ARM_DIR/$VCF_NAME.$REF.$ARM.snv.vcf.gz

    # Split by arm
    # Fill missing genotypes (with 0, i.e. reference allele, else msmc breaks)
    vcftools \
      --gzvcf $PATH_CHR_VCF_GZ \
      --bed <(echo -e "#chrom\tstart\tend"; grep -P "${REF}\t" ${ARM}s.bed) \
      --remove-indels --recode --keep-INFO-all --stdout |\
      python3 /lizardfs/guarracino/chromosome_communities/scripts/fill_genotypes.py diploid |
      #bcftools norm -f $PATH_REF_FASTA -c e -m - |\
      bgzip -@ 48 > $PATH_CHR_ARM_SNV_VCF_GZ
    tabix $PATH_CHR_ARM_SNV_VCF_GZ
    #cat out.log
  done
done
```

## MSMC

### Coalescence p-arms and q-arms separately and for each chromosome

Split variants by sample:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/msmc
cd /lizardfs/guarracino/chromosome_communities/msmc

(seq 1 22) | while read i; do
  REF=chr$i
  echo $REF 

  PATH_CHR_VCF_GZ=/lizardfs/erikg/HPRC/year1v2genbank/wgg.88/$REF.pan/$REF.pan.fa.a2fb268.4030258.6a1ecc2.smooth.chm13.vcf.gz
  VCF_NAME=$(basename $PATH_CHR_VCF_GZ .vcf.gz)

  for ARM in p_arm q_arm; do
    ARM_DIR=$REF/$ARM
    mkdir -p $ARM_DIR

    PATH_CHR_ARM_SNV_VCF_GZ=$ARM_DIR/$VCF_NAME.$REF.$ARM.snv.vcf.gz

    # Split by sample
    for SAMPLE in `bcftools query -l $PATH_CHR_ARM_SNV_VCF_GZ | tr '\n' ' '`; do
      echo $REF $ARM $SAMPLE
      
      # Fill missing genotypes (with 0, i.e. reference allele, else msmc breaks)
      PATH_CHR_ARM_SNV_SAMPLE_VCF_GZ=$ARM_DIR/$VCF_NAME.$REF.$ARM.$SAMPLE.snv.vcf.gz
      bcftools view --samples $SAMPLE $PATH_CHR_ARM_SNV_VCF_GZ |\
        bgzip -@ 48 -c > $PATH_CHR_ARM_SNV_SAMPLE_VCF_GZ
    done
  done
done
```


Collect a few samples (MSMC doesn't scale) to build the populations to compare:

```shell
NUM_SAMPLES=15

(seq 1 22) | while read i; do
  REF=chr$i
  echo $REF 

  PATH_CHR_VCF_GZ=/lizardfs/erikg/HPRC/year1v2genbank/wgg.88/$REF.pan/$REF.pan.fa.a2fb268.4030258.6a1ecc2.smooth.chm13.vcf.gz
  VCF_NAME=$(basename $PATH_CHR_VCF_GZ .vcf.gz)

  for ARM in p_arm q_arm; do
    ARM_DIR=$REF/$ARM
    PATH_CHR_ARM_SNV_VCF_GZ=$ARM_DIR/$VCF_NAME.$REF.$ARM.snv.vcf.gz
    PATH_CHR_ARM_POP_TXT=$ARM_DIR/$REF.$ARM.population.txt
    
    for SAMPLE in `bcftools query -l $PATH_CHR_ARM_SNV_VCF_GZ | grep grch38 -v | head -n ${NUM_SAMPLES} | tr '\n' ' '`; do
      echo $REF $ARM $SAMPLE
      
      PATH_CHR_ARM_SNV_SAMPLE_VCF_GZ=$ARM_DIR/$VCF_NAME.$REF.$ARM.$SAMPLE.snv.vcf.gz
      
      echo $PATH_CHR_ARM_SNV_SAMPLE_VCF_GZ >> $PATH_CHR_ARM_POP_TXT
    done
  done
done
```

Run `msmc`:

```shell
(seq 1 22) | while read i; do
  REF=chr$i
  echo $REF

  PATH_CHR_VCF_GZ=/lizardfs/erikg/.../wgg88/.../..$REF...vcf.gz
  VCF_NAME=$(basename $PATH_CHR_VCF_GZ .vcf.gz)

  for ARM in p_arm q_arm; do
    ARM_DIR=$REF/$ARM
    PATH_CHR_ARM_POP_TXT=$ARM_DIR/$REF.$ARM.population.txt

    sbatch -c 48 -p workers --job-name msmc-$i /lizardfs/guarracino/chromosome_communities/scripts/msmc2.sh $PATH_CHR_ARM_POP_TXT /home/guarracino/tools/msmc-tools/generate_multihetsep.py /home/guarracino/tools/msmc2_linux64bit $ARM_DIR/$REF.$ARM.msmc2
  done
done
```

### Coalescence p-arms together/all pairs/all combinations ???
...



## SMC++

### Coalescence p-arms and q-arms separately and for each chromosome

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/smcpp
cd /lizardfs/guarracino/chromosome_communities/smcpp

(seq 1 22) | while read i; do
  REF=chr$i
  echo $REF 

  for ARM in p_arm q_arm; do
    PATH_CHR_ARM_SNV_VCF_GZ=/lizardfs/guarracino/chromosome_communities/msmc/$REF/$ARM/$REF.pan.fa.a2fb268.4030258.6a1ecc2.smooth.chm13.$REF.$ARM.snv.vcf.gz
    SAMPLES=`bcftools query -l $PATH_CHR_ARM_SNV_VCF_GZ | grep grch38 -v | tr '\n' ',' | sed 's/.$//'`
    
    /home/guarracino/smcpp/bin/smc++ vcf2smc \
      $PATH_CHR_ARM_SNV_VCF_GZ \
      $REF.$ARM.smc.gz \
      chm13#$REF \
      $REF_$ARM:$SAMPLES
  done
done

(seq 1 22) | while read i; do
  REF=chr$i
  echo $REF 

  for ARM in p_arm q_arm; do
    PATH_CHR_ARM_SNV_VCF_GZ=/lizardfs/guarracino/chromosome_communities/msmc/$REF/$ARM/$REF.pan.fa.a2fb268.4030258.6a1ecc2.smooth.chm13.$REF.$ARM.snv.vcf.gz

    sbatch -c 48 -p workers --jon-name smcpp-$i --wrap'/home/guarracino/smcpp/bin/smc++ estimate 1e-8 '$REF'.'$ARM'.smc.gz -o '$REF'.'$ARM
  done
done

(seq 1 22) | while read i; do
  REF=chr$i
  echo $REF 

  #for ARM in p_arm q_arm; do
  #  PATH_CHR_ARM_SNV_VCF_GZ=/lizardfs/guarracino/chromosome_communities/msmc/$REF/$ARM/$REF.pan.fa.a2fb268.4030258.6a1ecc2.smooth.chm13.$REF.$ARM.snv.vcf.gz

  #  /home/guarracino/smcpp/bin/smc++ plot $REF.$ARM/$REF.$ARM.png $REF.$ARM/model.final.json
  #done
  /home/guarracino/smcpp/bin/smc++ plot $REF.pq_arms.png $REF.p_arm/model.final.json $REF.q_arm/model.final.json -g 1000

done





/home/guarracino/smcpp/bin/smc++ vcf2smc \
  /lizardfs/guarracino/chromosome_communities/msmc/chr21/p_arm/chr21.pan.fa.a2fb268.4030258.6a1ecc2.smooth.chm13.chr21.p_arm.snv.vcf.gz \
  chr21.smc.gz \
  chm13#chr21 \
  POP:HG00735,HG00741,HG01071,HG01106,HG01109,HG01123,HG01175,HG01243,HG01258,HG01358,HG01361,HG01891,HG01928,HG01952,HG01978,HG02055,HG02080,HG02109,HG02145,HG02148,HG02257,HG02486,HG02559,HG02572,HG02622



/home/guarracino/smcpp/bin/smc++ estimate 1e-8 chr21.smc.gz \
/home/guarracino/smcpp/bin/smc++ plot plot.png model.final.json  -g 1000
```


















Collect samples to build the populations to compare:

```shell
# Take pq-touching untangling samples
#(echo 13; echo 21) | while read i; do
#  REF=chr$i
#  
#  OUTPUT_DIR=13vs21_wrt_${REF}
#  mkdir -p $OUTPUT_DIR
#    
#  (echo 13; echo 21) | while read j; do
#    PATH_PQ_TOUCHING_TSV_GZ=/lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.untangle.chm13#chr$j.e50000.m1000.j08.n10.grounded.pq_touching.tsv.gz
#    SAMPLES_TO_TAKE=$(zcat $PATH_PQ_TOUCHING_TSV_GZ | sed '1d' | cut -f 1 | uniq | grep 'chm13\|grch38' -v)
#    
#    echo $SAMPLES_TO_TAKE | tr ' ' '\n' | head -n 2 | while read SAMPLE; do
#      echo "/lizardfs/guarracino/chromosome_communities/msmc/$REF/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.snv.${REF}.p_arm.${SAMPLE}.vcf.gz" \
#        >> $OUTPUT_DIR/$REF.samples.p_arms.txt
#      echo "/lizardfs/guarracino/chromosome_communities/msmc/$REF/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.snv.${REF}.q_arm.${SAMPLE}.vcf.gz" \
#        >> $OUTPUT_DIR/$REF.samples.q_arms.txt
#    done
#  done
#done

#All pq-mapping contigs is too much
#(seq 13 15; seq 21 22) | while read i; do
#  PATH_REF_FASTA_FAI=/lizardfs/guarracino/chromosome_communities/pq_contigs/chr$i.vs.chm13.50kbps.pq_contigs.union.hg002prox.fa.gz.fai
#
#  cat $PATH_REF_FASTA_FAI $nsamples | cut -f 1 | while read SAMPLE; do
#    echo "chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.snv.p_arms.${SAMPLE}.vcf.gz" >> samples.p_arms.txt
#    echo "chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.snv.q_arm.${SAMPLE}.vcf.gz" >> samples.q_arms.txt
#  done
#done

# Take pq-touching untangling samples
(seq 13 13) | while read i; do
  echo $i
  REF=chr$i
  
  OUTPUT_DIR=p${i}_vs_q${rm}_wrt_${REF}
  mkdir -p $OUTPUT_DIR
    
  PATH_PQ_TOUCHING_TSV_GZ=/lizardfs/guarracino/chromosome_communities/untangle/grounded/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.untangle.chm13#$REF.e50000.m1000.j08.n10.grounded.pq_touching.tsv.gz
  SAMPLES_TO_TAKE=$(zcat $PATH_PQ_TOUCHING_TSV_GZ | sed '1d' | cut -f 1 | uniq | grep 'chm13\|grch38' -v)
    
  echo $SAMPLES_TO_TAKE | tr ' ' '\n' | head -n 5 | while read SAMPLE; do
    echo "/lizardfs/guarracino/chromosome_communities/msmc/$REF/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.snv.${REF}.p_arm.${SAMPLE}.vcf.gz" \
      >> $OUTPUT_DIR/$REF.samples.p_arms.txt
    echo "/lizardfs/guarracino/chromosome_communities/msmc/$REF/chrACRO+refs.50kbps.pq_contigs.union.hg002prox.fa.gz.e998a33.4030258.fb5ffef.smooth.fix.chm13.snv.${REF}.q_arm.${SAMPLE}.vcf.gz" \
      >> $OUTPUT_DIR/$REF.samples.q_arms.txt
  done
done

(seq 13 13) | while read i; do
  echo $i
  REF=chr$i
  OUTPUT_DIR=p${i}_vs_q${rm}_wrt_${REF}
  
  bash /lizardfs/guarracino/chromosome_communities/scripts/msmc2pops.sh \
    $OUTPUT_DIR/$REF.samples.p_arms.txt \
    $OUTPUT_DIR/$REF.samples.q_arms.txt \
    /home/guarracino/tools/msmc-tools/generate_multihetsep.py \
    /home/guarracino/tools/msmc-tools/combineCrossCoal.py \
    /home/guarracino/tools/msmc2_linux64bit \
    /lizardfs/guarracino/chromosome_communities/msmc/$OUTPUT_DIR/$REF
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
```
