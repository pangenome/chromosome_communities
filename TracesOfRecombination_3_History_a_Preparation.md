# History

## Preparation

### Variables

```shell
THREADS=48
RUN_MINIMAP2=/gnu/store/7xpn4bkg6jknk3mzdk0alkrxz5i40j8c-minimap2-2.24/bin/minimap2

# HACKED bamCaller.py to emit haploid VCF files --> if altAllele != "." and re.match("^[01]", genotypes):
RUN_HACKED_BAM_CALLER=/home/guarracino/tools/msmc-tools/bamCaller.py

# Using a hacked version, changing a few rows in generate_multihetsep.py to get phased fake-diploid data
#if len(geno) != 3 or geno[1] not in ['|/']:
#  print ("HACKED: Non-diploid SNP found and considered as phased data: %s" % geno, file=sys.stderr)
#  phased = True
#  geno = "%s|%s" % (geno[0], geno[0])
RUN_HACKED_GENERATE_MULTIHETSEP=/home/guarracino/tools/msmc-tools/generate_multihetsep.py

RUN_COMBINE_CROSS_COAL=/home/guarracino/tools/msmc-tools/combineCrossCoal.py
RUN_MSMC2=/home/guarracino/tools/msmc2_linux64bit
```

### Steps

Split `chrS` by haplotype:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/haplotypes
cd /lizardfs/guarracino/chromosome_communities/haplotypes

cat /lizardfs/erikg/HPRC/year1v2genbank/parts/chrS.pan.fa.fai | cut -f 1,2 -d '#' | grep chr -v | sort | uniq > haplotypes.txt

( seq 13 15; seq 21 22 ) | while read i; do
  mkdir -p chr$i
  
  PATH_CHR_FASTA=/lizardfs/erikg/HPRC/year1v2genbank/parts/chr$i.pan.fa
  cat haplotypes.txt | while read HAPLO; do
    echo "chr$i" $HAPLO
  
    samtools faidx $PATH_CHR_FASTA $(grep "^$HAPLO" $PATH_CHR_FASTA.fai | cut -f 1) |
      bgzip -@ THREADS -c > chr$i/chr$i.$HAPLO.fa.gz
    samtools faidx chr$i/chr$i.$HAPLO.fa.gz
  done
  
  PATH_CHR_FASTA=/lizardfs/guarracino/chromosome_communities/assemblies/hg002.chr$i.prox.fa.gz
  for HAPLO in HG002#MAT HG002#PAT; do
    echo "chr$i" $HAPLO
    
    samtools faidx $PATH_CHR_FASTA $(grep "^$HAPLO" $PATH_CHR_FASTA.fai | cut -f 1) |
      bgzip -@ THREADS -c > chr$i/chr$i.$HAPLO.fa.gz
    samtools faidx chr$i/chr$i.$HAPLO.fa.gz
  done
done
```

Map contigs against acrocentric chromosomes:

```shell
( seq 13 15; seq 21 22 ) | while read i; do    
  ( seq 13 15; seq 21 22 ) | while read j; do
    PATH_REF_CHR_FASTA=/lizardfs/guarracino/chromosome_communities/assemblies/chm13.chr$j.fa.gz

    ( cat haplotypes.txt; echo HG002#MAT; echo HG002#PAT ) | while read HAPLO; do
      echo "chr$i" $HAPLO
      #$RUN_MINIMAP2 $PATH_REF_CHR_FASTA chr$i/chr$i.$HAPLO.fa.gz -t $THREADS -x asm20 -R "@RG\tID:$HAPLO\tSM:$HAPLO" -ca | samtools sort > chr$i/chr$i.$HAPLO.on.chr$j.bam
      sbatch -p workers -c 2 --job-name "$i-on-$j" --wrap 'hostname; \time -v '$RUN_MINIMAP2' '$PATH_REF_CHR_FASTA' /lizardfs/guarracino/chromosome_communities/haplotypes/chr'$i'/chr'$i'.'$HAPLO'.fa.gz -t 2 -x asm20 -R "@RG\tID:'$HAPLO'\tSM:'$HAPLO'" -ca | samtools sort > /lizardfs/guarracino/chromosome_communities/haplotypes/chr'$i'/chr'$i'.'$HAPLO'.on.chr'$j'.bam'
    done
  done
done
```

Generating VCF and Mask files:

```shell
( seq 13 15; seq 21 22 ) | while read i; do  
  ( seq 13 15; seq 21 22 ) | while read j; do
    PATH_REF_CHR_FASTA=/lizardfs/guarracino/chromosome_communities/assemblies/chm13.chr$j.fa.gz
      
    ( cat haplotypes.txt; echo HG002#MAT; echo HG002#PAT ) | while read HAPLO; do
      echo "chr$i" $HAPLO
      PREFIX=/lizardfs/guarracino/chromosome_communities/haplotypes/chr$i/chr$i.$HAPLO.on.chr$j
      
      #bcftools mpileup -B -q 0 -Q 0 -C 50 -f $PATH_REF_CHR_FASTA $PREFIX.bam |\
      #  bcftools call -c --skip-variants indels --ploidy 1 --threads $THREADS |\
      #  $RUN_HACKED_BAM_CALLER 1 $PREFIX.mask.bed.gz --minMapQ 1 |\
      #  bgzip -@ $THREADS -c > $PREFIX.vcf.gz
      sbatch -p workers -c 2 --job-name "$i-on-$j" --wrap 'hostname; \time -v bcftools mpileup -B -q 0 -Q 0 -C 50 -f '$PATH_REF_CHR_FASTA' '$PREFIX'.bam | bcftools call -c --skip-variants indels --ploidy 1 --threads 2 | '$RUN_HACKED_BAM_CALLER' 1 '$PREFIX'.mask.bed.gz --minMapQ 1 | bgzip -@ 2 -c > '$PREFIX'.vcf.gz'
    done
  done
done
```

Split variants by arm:

```shell
# Prepare p/q-arms BED files:
sed 's/chr/chm13#chr/g' /lizardfs/guarracino/chromosome_communities/data/chm13.centromeres.approximate.bed | \
  bedtools sort | \
  bedtools complement \
    -i - \
    -g <(cut -f 1,2 /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa.fai | grep 'chrM\|chrX' -v | sort) \
  > tmp.bed
  
# Take odd rows
(echo -e "#chrom\tstart\tend"; sed -n 1~2p tmp.bed) > p_arms.bed
# Take even rows
(echo -e "#chrom\tstart\tend"; sed -n 2~2p tmp.bed) > q_arms.bed
rm tmp.bed

# Split
( seq 13 15; seq 21 22 ) | while read i; do  
  ( seq 13 15; seq 21 22 ) | while read j; do
    PATH_REF_CHR_FASTA=/lizardfs/guarracino/chromosome_communities/assemblies/chm13.chr$j.fa.gz
      
    ( cat haplotypes.txt; echo HG002#MAT; echo HG002#PAT ) | while read HAPLO; do
      echo "chr$i" $HAPLO
      PREFIX=/lizardfs/guarracino/chromosome_communities/haplotypes/chr$i/chr$i.$HAPLO.on.chr$j
      
      PATH_VCF_GZ=$PREFIX.vcf.gz
      PATH_MASK_BED_GZ=$PREFIX.mask.bed.gz
      REF=chm13#chr$j
      
      for ARM in p_arm q_arm; do
        PATH_ARM_VCF_GZ=$PREFIX.$ARM.vcf.gz
        vcftools \
          --gzvcf $PATH_VCF_GZ \
          --bed <(echo -e "#chrom\tstart\tend"; grep -P "${REF}\t" ${ARM}s.bed) \
          --recode --keep-INFO-all --out chr$i.$HAPLO.on.chr$j.$ARM --stdout |\
          #bcftools norm -f $PATH_REF_CHR_FASTA -c e -m - |\
          bgzip -@ $THREADS > $PATH_ARM_VCF_GZ
        tabix --force $PATH_ARM_VCF_GZ
        
        PATH_ARM_MASK_BED_GZ=$PREFIX.$ARM.mask.bed.gz
        bedtools intersect \
          -a <(zcat $PATH_MASK_BED_GZ) \
          -b <(grep -P "${REF}\t" ${ARM}s.bed) |\
          bgzip -@ $THREADS > $PATH_ARM_MASK_BED_GZ
      done
    done
  done
done
```

Haplotypes to ignore for the p-arms:

```shell
# Consider only variants on the right of the rDNA
zgrep rDNA /lizardfs/guarracino/chromosome_communities/data/chm13.CentromeresAndTelomeres.CenSatAnnotation.txt.gz | \
  sed 's/chr/chm13#chr/g' | \
  awk -v OFS='\t' '{print $1,"0",$3}' |\
  bedtools complement \
    -i - \
    -g <(grep 'chr13\|chr14\|chr15\|chr21\|chr22' /lizardfs/erikg/HPRC/year1v2genbank/assemblies/chm13.fa.fai | cut -f 1,2 | sort) \
  > rDNA.right.bed
  
( seq 13 15; seq 21 22 ) | while read i; do  
  ( seq 13 15; seq 21 22 ) | while read j; do
    PATH_REF_CHR_FASTA=/lizardfs/guarracino/chromosome_communities/assemblies/chm13.chr$j.fa.gz
      
    ( cat haplotypes.txt; echo HG002#MAT; echo HG002#PAT ) | while read HAPLO; do
      echo "chr$i" $HAPLO
      PREFIX=/lizardfs/guarracino/chromosome_communities/haplotypes/chr$i/chr$i.$HAPLO.on.chr$j
      
      PATH_ARM_VCF_GZ=$PREFIX.p_arm.vcf.gz
      PATH_ARM_MASK_BED_GZ=$PREFIX.p_arm.mask.bed.gz
      REF=chm13#chr$j
      
      PATH_ARM_RIGHT_VCF_GZ=$PREFIX.p_arm.right.vcf.gz
      vcftools \
        --gzvcf $PATH_ARM_VCF_GZ \
        --bed <(echo -e "#chrom\tstart\tend"; grep -P "${REF}\t" rDNA.right.bed) \
        --recode --keep-INFO-all --out chr$i.$HAPLO.on.chr$j.p_arm.right --stdout |\
        #bcftools norm -f $PATH_REF_CHR_FASTA -c e -m - |\
        bgzip -@ $THREADS > $PATH_ARM_RIGHT_VCF_GZ
      tabix --force $PATH_ARM_RIGHT_VCF_GZ
      
      PATH_ARM_RIGHT_MASK_BED_GZ=$PREFIX.p_arm.right.mask.bed.gz
      bedtools intersect \
        -a <(zcat $PATH_ARM_MASK_BED_GZ) \
        -b <(grep -P "${REF}\t" rDNA.right.bed) |\
        bgzip -@ $THREADS > $PATH_ARM_RIGHT_MASK_BED_GZ
    done
  done
done
```

Combining into one file multiple individuals from all possible pairs of populations:

```shell
mkdir -p  /lizardfs/guarracino/chromosome_communities/haplotypes/msmc

#( seq 13 15; seq 21 22 ) | while read r; do
( echo 13 ) | while read r; do
  # Skip haplotypes with less than 5000 called variants on the p-arm
  grep 'Sites' *on.chr$r.p_arm.right.log | cut -f 1,4 -d ' ' | awk '{if ($2 < 5000) print($0)}' | sort -k 2nr | cut -f 2 -d '.' | sort | uniq > haplotypes_to_skip.on.chr$r.p_arm.right.txt

  # All possible acro-pairs
  ( seq 13 15; seq 21 22 ) | while read a; do 
    ( seq 13 15; seq 21 22 ) | while read b; do
      if (( $a < $b )); then
        echo "chr$a vs chr$b on chr$r"
        
        #MASKsA=$( ls /lizardfs/guarracino/chromosome_communities/haplotypes/chr$a/chr$a.*.on.chr$r.p_arm.right.mask.bed.gz | grep -f haplotypes_to_skip.on.chr$r.p_arm.right.txt -v | sed -z 's/\n/ --mask /g' | sed 's/ --mask $//g')
        #VCFsA=$( ls /lizardfs/guarracino/chromosome_communities/haplotypes/chr$a/chr$a.*.on.chr$r.p_arm.right.vcf.gz | grep -f haplotypes_to_skip.on.chr$r.p_arm.right.txt -v | tr '\n' ' ')
        #MASKsB=$( ls /lizardfs/guarracino/chromosome_communities/haplotypes/chr$b/chr$b.*.on.chr$r.p_arm.right.mask.bed.gz | grep -f haplotypes_to_skip.on.chr$r.p_arm.right.txt -v | sed -z 's/\n/ --mask /g' | sed 's/ --mask $//g')
        #VCFsB=$( ls /lizardfs/guarracino/chromosome_communities/haplotypes/chr$b/chr$b.*.on.chr$r.p_arm.right.vcf.gz | grep -f haplotypes_to_skip.on.chr$r.p_arm.right.txt -v | tr '\n' ' ')
        #sbatch -p workers -c 2 --job-name "$a-$b-$r" --wrap 'hostname; \time -v '$RUN_HACKED_GENERATE_MULTIHETSEP' --mask '$MASKsA' --mask '$MASKsB' '$VCFsA' '$VCFsB' > /lizardfs/guarracino/chromosome_communities/haplotypes/multihetsep/chr'$a'.vs.chr'$b'.on.chr'$r'.phased.hetsep'
        
        PATH_POPULATION_1_TXT=/lizardfs/guarracino/chromosome_communities/haplotypes/msmc/chr$a.on.chr$r.txt
        ls /lizardfs/guarracino/chromosome_communities/haplotypes/chr$a/chr$a.*.on.chr$r.p_arm.right.vcf.gz | grep -f haplotypes_to_skip.on.chr$r.p_arm.right.txt -v > $PATH_POPULATION_1_TXT
        PATH_POPULATION_2_TXT=/lizardfs/guarracino/chromosome_communities/haplotypes/msmc/chr$b.on.chr$r.txt
        ls /lizardfs/guarracino/chromosome_communities/haplotypes/chr$b/chr$b.*.on.chr$r.p_arm.right.vcf.gz | grep -f haplotypes_to_skip.on.chr$r.p_arm.right.txt -v > $PATH_POPULATION_2_TXT
        
        sbatch -c 48 -p allhighmem --job-name $a-$b-$r /lizardfs/guarracino/chromosome_communities/scripts/msmc2pops.sh $PATH_POPULATION_1_TXT $PATH_POPULATION_2_TXT $RUN_HACKED_GENERATE_MULTIHETSEP $RUN_COMBINE_CROSS_COAL $RUN_MSMC2 chr$a.vs.chr$b.on.chr$r
      fi
    done
  done
done


```