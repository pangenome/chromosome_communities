## Recombination rate

```shell
###################################
# Octopus
THREADS=48
pathLDhat=~/tools/LDhat
runPhi=~/tools/PhiPack/src/Phi
runSumsLDhat=/lizardfs/guarracino/chromosome_communities/scripts/Sums_LDHat_pack.sh
pathLDhatChunk=/lizardfs/guarracino/chromosome_communities/scripts/LDhatChunk.sh
pathGetIndexesR=/lizardfs/guarracino/chromosome_communities/scripts/get_indexes.R
pathGetRecombinationRatePlotR=/lizardfs/guarracino/chromosome_communities/scripts/get_recombination_rate_plot.R

OUTPUT_DIR=/lizardfs/guarracino/chromosome_communities/recombination_rate
####################

###################################
# Locally (for testing)
THREADS=16
pathLDhat=~/git/LDhat/
runPhi=~/Downloads/Pangenomics/LDhat/PhiPack/src/Phi
runSumsLDhat=~/git/chromosome_communities/scripts/Sums_LDHat_pack.sh
pathLDhatChunk=~/git/chromosome_communities/scripts/LDhatChunk.sh
pathGetIndexesR=~/git/chromosome_communities/scripts/get_indexes.R
pathGetRecombinationRatePlotR=~/git/chromosome_communities/scripts/get_recombination_rate_plot.R

OUTPUT_DIR=~/Downloads/Pangenomics/LDhat/recombination_rate
###################################
```

Call variants in a haploid setting:

```shell
sbatch -p workers -c 48 --wrap 'vg deconstruct -P chm13 -H '?' -e -a -t 48 /lizardfs/guarracino/chromosome_communities/graphs/chrA.pan.s100k.l300k.p97.n200/chrA.pan.fa.4825454.4030258.0517e7e.smooth.fix.gfa | bgzip -@ 48 > /lizardfs/guarracino/chromosome_communities/graphs/chrA.pan.s100k.l300k.p97.n200/chrA.pan.fa.4825454.4030258.0517e7e.smooth.fix.chm13.vcf.gz && tabix /lizardfs/guarracino/chromosome_communities/graphs/chrA.pan.s100k.l300k.p97.n200/chrA.pan.fa.4825454.4030258.0517e7e.smooth.fix.chm13.vcf.gz'
sbatch -p workers -c 48 --wrap 'vg deconstruct -P grch38 -H '?' -e -a -t 48 /lizardfs/guarracino/chromosome_communities/graphs/chrA.pan.s100k.l300k.p97.n200/chrA.pan.fa.4825454.4030258.0517e7e.smooth.fix.gfa | bgzip -@ 48 > /lizardfs/guarracino/chromosome_communities/graphs/chrA.pan.s100k.l300k.p97.n200/chrA.pan.fa.4825454.4030258.0517e7e.smooth.fix.grch38.vcf.gz && tabix /lizardfs/guarracino/chromosome_communities/graphs/chrA.pan.s100k.l300k.p97.n200/chrA.pan.fa.4825454.4030258.0517e7e.smooth.fix.grch38.vcf.gz'

sbatch -p workers -c 48 --wrap 'vg deconstruct -P chm13 -H '?' -e -a -t 48 /lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.100kbps.pq_contigs.union.s100k.l300k.p98.n97/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.gfa | bgzip -@ 48 > /lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.100kbps.pq_contigs.union.s100k.l300k.p98.n97/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.chm13.vcf.gz && tabix /lizardfs/guarracino/chromosome_communities/graphs/chrA.pan.s100k.l300k.p97.n200/chrA.pan.fa.4825454.4030258.0517e7e.smooth.fix.chm13.vcf.gz'
sbatch -p workers -c 48 --wrap 'vg deconstruct -P grch38 -H '?' -e -a -t 48 /lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.100kbps.pq_contigs.union.s100k.l300k.p98.n97/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.gfa | bgzip -@ 48 > /lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.100kbps.pq_contigs.union.s100k.l300k.p98.n97/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.grch38.vcf.gz && tabix /lizardfs/guarracino/chromosome_communities/graphs/chrA.pan.s100k.l300k.p97.n200/chrA.pan.fa.4825454.4030258.0517e7e.smooth.fix.grch38.vcf.gz'
```

Prepare the VCF file with only SNPs:

```shell
PATH_VCF=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.100kbps.pq_contigs.union.s100k.l300k.p98.n97/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.chm13.vcf.gz
#PATH_VCF=~/Downloads/Pangenomics/LDhat/SimulatedPopulations/HatLandscapeN16Len1000000Nrhs15_th0.01_540_1.pggb/HatLandscapeN16Len1000000Nrhs15_th0.01_540_1.fixed.fa.5c3c9a3.7bdde5a.a933754.smooth.fix.gfa.vcf.gz

VCF_NAME=$(basename $PATH_VCF .vcf.gz)
SEG_LENGTH=1000 # 1kb is recommended for LDJump


mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

# Get reference names
zgrep '^#' $PATH_VCF -v | cut -f 1 | sort | uniq > ref_names.txt

# guix install vcftools
cat ref_names.txt | while read REF_NAME; do
  echo $REF_NAME
  
  #PATH_REF_FASTA=/lizardfs/guarracino/chromosome_communities/assemblies/$(echo $REF_NAME | tr '#' '.').fa.gz
  PATH_REF_FASTA=~/Downloads/Pangenomics/LDhat/SimulatedPopulations/HatLandscapeN16Len1000000Nrhs15_th0.01_540_1.fixed.$REF_NAME.fa

  PATH_VCF_SNPS=$OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.vcf.gz

  # Take SNPs and normalize
  vcftools --gzvcf $PATH_VCF --chr ${REF_NAME} --remove-indels --recode --recode-INFO-all --stdout |\
    bcftools norm -f $PATH_REF_FASTA -c s -m - |\
    bgzip -@ $THREADS > $OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.tmp.vcf.gz
  rm out.log
  
  # Take only samples with at least 1 variants in the VCF file
  sample_with_variants=$(bcftools stats -s '-' $OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.tmp.vcf.gz | awk '$1=="PSC" && $12+$13>0 {print $3}' | paste -s -d',')
  bcftools view --samples $sample_with_variants $OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.tmp.vcf.gz |\
   bgzip -@ $THREADS > $PATH_VCF_SNPS
  tabix $PATH_VCF_SNPS

  rm $OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.tmp.vcf.gz
done
```

In parallel:
- prepare range;
- chunk VCF file with the SNPs
- get a FASTA by applying for each sample the variants to the reference
- apply LDhat on the FASTA
- get indexes for the FASTA

```shell
cat ref_names.txt | while read REF_NAME; do
    echo $REF_NAME
    
    PATH_REF_FASTA=/lizardfs/guarracino/chromosome_communities/assemblies/$(echo $REF_NAME | tr '#' '.').fa.gz
    #PATH_REF_FASTA=~/Downloads/Pangenomics/LDhat/SimulatedPopulations/HatLandscapeN16Len1000000Nrhs15_th0.01_540_1.fixed.$REF_NAME.fa
    
    LENGTH_OF_SEQ=$(cut -f 2 ${PATH_REF_FASTA}.fai)
    START_OF_SEQ=0
    END_OF_SEQ=$LENGTH_OF_SEQ
    
    NUM_OF_SEGS=$(echo $LENGTH_OF_SEQ / $SEG_LENGTH | bc)
    NUM_OF_SEGS_PLUS_1=$(echo $LENGTH_OF_SEQ / $SEG_LENGTH + 1 | bc)
    
    PATH_VCF_SNPS=$OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.vcf.gz
    VCF_SNPS_NAME=$(basename $PATH_VCF_SNPS .vcf.gz)
    
    N=$(zgrep '^#CHROM' $PATH_VCF_SNPS -m 1 | cut -f 10- | tr '\t' '\n' | grep grch38 -v | wc -l)
    
    MIN_POS=$(zgrep '^#' $PATH_VCF_SNPS -v | cut -f 2 | sort -k 1n | head -n 1)
    MAX_POS=$(zgrep '^#' $PATH_VCF_SNPS -v | cut -f 2 | sort -k 1n | tail -n 1)

    cd /scratch/
    mkdir -p ${REF_NAME}-recomb
    TEMP_DIR=${REF_NAME}-recomb/$REF_NAME-temp
    mkdir -p $TEMP_DIR

    seq 1 $NUM_OF_SEGS_PLUS_1 | parallel -j $THREADS bash $pathLDhatChunk {} \
      $SEG_LENGTH \
      $NUM_OF_SEGS_PLUS_1 \
      $END_OF_SEQ \
      $N \
      $MIN_POS \
      $MAX_POS \
      $PATH_REF_FASTA \
      $REF_NAME \
      $TEMP_DIR \
      $PATH_VCF_SNPS \
      $runSumsLDhat \
      $pathLDhat \
      $runPhi \
      $pathGetIndexesR
    
    INDEXES_DIR=$OUTPUT_DIR/$REF_NAME-indexes
    mkdir -p $INDEXES_DIR
    
    # With cat we can have 'Argument list too long'
    # Retrieve files in order
    for file in $(find ${TEMP_DIR} -name '*.indexes.*.tsv' -type f | sort); do
        cat $file >> ${INDEXES_DIR}/indexes.tsv
    done
    for file in $(find ${TEMP_DIR} -name 'Sums_part_main_job*.txt' -type f | sort); do
        cat $file >> ${INDEXES_DIR}/Sums_part_main_job.txt
    done
    
    # Cleaning the working directory
    rm Phi*
    rm out.log
done
```

Recombination plots:

```shell
cat ref_names.txt | while read REF_NAME; do
    echo $REF_NAME
    
    PATH_REF_FASTA=/lizardfs/guarracino/chromosome_communities/assemblies/$(echo $REF_NAME | tr '#' '.').fa.gz
    #PATH_REF_FASTA=~/Downloads/Pangenomics/LDhat/SimulatedPopulations/HatLandscapeN16Len1000000Nrhs15_th0.01_540_1.fixed.$REF_NAME.fa
    
    LENGTH_OF_SEQ=$(cut -f 2 ${PATH_REF_FASTA}.fai)
    
    PATH_VCF_SNPS=$OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.vcf.gz
    N=$(zgrep '^#CHROM' $PATH_VCF_SNPS -m 1 | cut -f 10- | tr '\t' '\n' | grep grch38 -v | wc -l)
    
    INDEXES_DIR=$OUTPUT_DIR/$REF_NAME-indexes

    # guix install r-data-table
    Rscript $pathGetRecombinationRatePlotR ${INDEXES_DIR} $SEG_LENGTH $N $LENGTH_OF_SEQ $REF_NAME
done
```

```shell

####cat $VCF|  vcflib vcfsnps | bgzip -c > $VCF_NAME.snps.vcf.gz && tabix $VCF_NAME.snps.vcf.gz
#vcftools --gzvcf $VCF --remove-indels --recode --recode-INFO-all --stdout | bcftools norm -f chm13.ACRO.fa -c s -m - | bgzip -@ 48 > $VCF_NAME.snps.vcf.gz
#tabix $VCF_NAME.snps.vcf.gz
##Missing header: zgrep '^#' -v $VCF_NAME.snps.vcf.gz | cut -f 1 | uniq | xargs -i echo tabix $VCF_NAME.snps.vcf.gz {} \| bgzip \> $VCF_NAME.snps.{}.vcf.gz \&|sh

# guix install vcftools
for i in 13 14 15 21 22; do
  chr=chm13#chr$i
  echo $chr

  vcftools --gzvcf $VCF --chr $chr --remove-indels --recode --recode-INFO-all --stdout |\
    bcftools norm -f /lizardfs/guarracino/chromosome_communities/assemblies/chm13.chr$i.fa.gz -c s -m - |\
    bgzip -@ 48 > $VCF_NAME.$chr.snps.tmp.vcf.gz

  sample_with_variants=$(bcftools stats -s '-' $VCF_NAME.$chr.snps.tmp.vcf.gz | awk '$1=="PSC" && $12+$13>0 {print $3}' | paste -s -d',')
  bcftools view --samples $sample_with_variants $VCF_NAME.$chr.snps.tmp.vcf.gz | bgzip -@ 48 > $VCF_NAME.$chr.snps.vcf.gz
  tabix $VCF_NAME.$chr.snps.vcf.gz
  
  rm $VCF_NAME.$chr.snps.tmp.vcf.gz
done

for i in 13 14 15 21 22; do
  chr=chm13#chr$i
  echo $chr
  
  zgrep '^#CHROM' $VCF_NAME.$chr.snps.vcf.gz -m 1 | cut -f 10- | tr '\t' '\n' | grep grch38 -v | while read SAMPLE; do
    echo $SAMPLE
    bcftools consensus -s $SAMPLE -f /lizardfs/guarracino/chromosome_communities/assemblies/chm13.chr$i.fa.gz $VCF_NAME.$chr.snps.vcf.gz |\
      sed "s/$chr/$SAMPLE/g" >> $VCF_NAME.$chr.snps.fa  
  done
  
  bgzip -@ 48 $VCF_NAME.$chr.snps.fa
  samtools faidx $VCF_NAME.$chr.snps.fa.gz
done


#guix install r-remotes
#guix install dos2unix

```

