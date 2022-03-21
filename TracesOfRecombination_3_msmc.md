# msmc

## Tools

```shell
cd ~/tools
git clone --recursive https://github.com/stschiff/msmc-tools.git

wget -c https://github.com/stschiff/msmc/releases/download/v1.1.0/msmc_1.1.0_linux64bit
chmod +x msmc_1.1.0_linux64bit
```

## Steps

Tutorial: https://github.com/stschiff/msmc-tools/blob/master/msmc-tutorial/guide.md#estimating-population-separation-history

Prepare VCF files for each sample, with respect to each reference:

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
  
# Take odd rows
(echo -e "#chrom\tstart\tend"; sed -n 1~2p tmp.bed) > p_arms.bed
# Take even rows
(echo -e "#chrom\tstart\tend"; sed -n 2~2p tmp.bed) > q_arms.bed
rm tmp.bed

(seq 13 15; seq 21 22) | while read i; do
  REF=chr$i
  echo $REF
    
  mkdir -p $REF
    
  # NOTE: in vctools, the BED file is expected to have a header line.
  vcftools \
    --gzvcf $PATH_SNV_VCF_GZ \
    --bed <(echo -e "#chrom\tstart\tend"; grep $REF p_arms.bed) \
    --recode --out $REF.p --keep-INFO-all --stdout |\
    #bcftools norm -f $PATH_REF_FASTA -c e -m - |\
    bgzip -@ $THREADS > $REF/$VCF_NAME.$REF.p_arm.vcf.gz
  vcftools \
    --gzvcf $PATH_SNV_VCF_GZ \
    --bed <(echo -e "#chrom\tstart\tend"; grep $REF q_arms.bed) \
    --recode --out $REF.q --keep-INFO-all --stdout |\
    #bcftools norm -f $PATH_REF_FASTA -c e -m - |\
    bgzip -@ $THREADS > $REF/$VCF_NAME.$REF.q_arm.vcf.gz
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
        bgzip -@ $THREADS -c > $REF/$VCF_NAME.$REF.$ARM.$SAMPLE.vcf.gz
    done
  done
done
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
  
  bash /lizardfs/guarracino/chromosome_communities/scripts/msmc.sh \
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