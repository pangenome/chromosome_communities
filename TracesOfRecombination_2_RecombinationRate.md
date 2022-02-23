## Recombination rate

```shell
# Octopus#################################################
THREADS=48
pathLDhat=~/tools/LDhat
runPhi=~/tools/PhiPack/src/Phi
runSumsLDhat=/lizardfs/guarracino/chromosome_communities/scripts/Sums_LDHat_pack.sh
pathLDhatChunk=/lizardfs/guarracino/chromosome_communities/scripts/LDhatChunk.sh
pathGetIndexesR=/lizardfs/guarracino/chromosome_communities/scripts/get_indexes.R
pathGetRecombinationRatePlotR=/lizardfs/guarracino/chromosome_communities/scripts/get_recombination_rate_plot.R

OUTPUT_DIR=/lizardfs/guarracino/chromosome_communities/recombination_rate
###################################

# Locally (for testing) ##################################
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
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

PATH_VCF=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.100kbps.pq_contigs.union.s100k.l300k.p98.n97/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.chm13.vcf.gz
#PATH_VCF=~/Downloads/Pangenomics/LDhat/SimulatedPopulations/HatLandscapeN16Len1000000Nrhs15_th0.01_540_1.pggb/HatLandscapeN16Len1000000Nrhs15_th0.01_540_1.fixed.fa.5c3c9a3.7bdde5a.a933754.smooth.fix.gfa.vcf.gz
#PATH_VCF=/lizardfs/guarracino/chromosome_communities/HatLandscapeN16Len1000000Nrhs15_th0.01_540_1.fixed.fa.5c3c9a3.7bdde5a.a933754.smooth.fix.gfa.vcf.gz

VCF_NAME=$(basename $PATH_VCF .vcf.gz)
SEG_LENGTH=1000 # 1kb is recommended for LDJump


# Get reference names
zgrep '^#' $PATH_VCF -v | cut -f 1 | sort | uniq > ref_names.txt

# guix install vcftools
cat ref_names.txt | while read REF_NAME; do
  echo $REF_NAME
  
  PATH_REF_FASTA=/lizardfs/guarracino/chromosome_communities/assemblies/$(echo $REF_NAME | tr '#' '.').fa.gz

  PATH_VCF_SNPS=$OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.vcf.gz

  # Take SNPs and normalize
  vcftools --gzvcf $PATH_VCF --chr ${REF_NAME} --remove-indels --recode --recode-INFO-all --stdout |\
    bcftools norm -f $PATH_REF_FASTA -c s -m - |\
    bgzip -@ $THREADS > $OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.tmp.vcf.gz
  rm out.log
  
  # Take only samples with at least 1 genotype in the VCF file
  sample_with_variants=$(bcftools stats -s '-' $OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.tmp.vcf.gz | awk '$1=="PSC" && $12+$13>0 {print $3}' | paste -s -d',')
  bcftools view --samples $sample_with_genotypes $OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.tmp.vcf.gz |\
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

    LENGTH_OF_SEQ=$(cut -f 2 ${PATH_REF_FASTA}.fai)
    START_OF_SEQ=0
    END_OF_SEQ=$LENGTH_OF_SEQ
    
    NUM_OF_SEGS=$(echo $LENGTH_OF_SEQ / $SEG_LENGTH | bc)
    NUM_OF_SEGS_PLUS_1=$(echo $LENGTH_OF_SEQ / $SEG_LENGTH + 1 | bc)
    
    PATH_VCF_SNPS=$OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.vcf.gz
    
    N=$(zgrep '^#CHROM' $PATH_VCF_SNPS -m 1 | cut -f 10- | tr '\t' '\n' | grep grch38 -v | wc -l)
 
    #Unused   
    #VCF_SNPS_NAME=$(basename $PATH_VCF_SNPS .vcf.gz)
    #zgrep '^#' $PATH_VCF_SNPS -v | cut -f 2 | sort -k 1n > $VCF_SNPS_NAME.all.positions.txt
   
    cd /scratch/
    TEMP_DIR=${REF_NAME}_s${SEG_LENGTH}-temp
    mkdir -p ${REF_NAME}_s${SEG_LENGTH}-recomb/${TEMP_DIR}
    cd ${REF_NAME}_s${SEG_LENGTH}-recomb

    seq 1 $NUM_OF_SEGS_PLUS_1 | parallel -j $THREADS bash $pathLDhatChunk {} \
      $SEG_LENGTH \
      $NUM_OF_SEGS_PLUS_1 \
      $END_OF_SEQ \
      $N \
      $PATH_REF_FASTA \
      $REF_NAME \
      $TEMP_DIR \
      $PATH_VCF_SNPS \
      $runSumsLDhat \
      $pathLDhat \
      $runPhi \
      $pathGetIndexesR
done
```

Retrieve chunk outputs and put all together in two matrices plus a file with additional information:

```shell
cat ref_names.txt | while read REF_NAME; do
    echo $REF_NAME

    cd /scratch/
    TEMP_DIR=${REF_NAME}_s${SEG_LENGTH}-recomb/${REF_NAME}_s${SEG_LENGTH}-temp
    
    INDEXES_DIR=$OUTPUT_DIR/${REF_NAME}_s${SEG_LENGTH}-indexes
    mkdir -p $INDEXES_DIR
    
    # With cat we can have 'Argument list too long'
    # Retrieve files in order
    for file in $(find ${TEMP_DIR} -name "${REF_NAME}.indexes.*.tsv" -type f | sort); do
        cat $file >> Indexes.tsv
    done
    mv Indexes.tsv ${INDEXES_DIR}/
    
    for file in $(find ${TEMP_DIR} -name "$REF_NAME.sums_part_main.*.txt" -type f | sort); do
        cat $file >> Sums_part_main_job.txt
    done
    mv Sums_part_main_job.txt ${INDEXES_DIR}/
    
    PATH_REF_FASTA=/lizardfs/guarracino/chromosome_communities/assemblies/$(echo $REF_NAME | tr '#' '.').fa.gz
    PATH_VCF_SNPS=$OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.vcf.gz
    N=$(zgrep '^#CHROM' $PATH_VCF_SNPS -m 1 | cut -f 10- | tr '\t' '\n' | grep grch38 -v | wc -l)
    LENGTH_OF_SEQ=$(cut -f 2 ${PATH_REF_FASTA}.fai)
    printf $SEG_LENGTH"\t"$N"\t"$LENGTH_OF_SEQ"\n" > ${INDEXES_DIR}/Information.txt
done
```

Recombination plots:

```shell
cat ref_names.txt | while read REF_NAME; do
    echo $REF_NAME
        
    #INDEXES_DIR=$OUTPUT_DIR/${REF_NAME}_s${SEG_LENGTH}-indexes
    INDEXES_DIR=~/Downloads/Pangenomics/LDhat/recombination_rate/${REF_NAME}_s${SEG_LENGTH}-indexes

    # guix install r-data-table
    Rscript $pathGetRecombinationRatePlotR \
      ${INDEXES_DIR}/Indexes.tsv \
      ${INDEXES_DIR}/Sums_part_main_job.txt \
      ${INDEXES_DIR}/Information.txt \
      ${INDEXES_DIR}/$REF_NAME \
      $REF_NAME
done
```

### Bootstrap

Sample 1 contig for each acrocentric chromosomes, for 1000 times:

```shell
# Sample 1 contig for each acros
NUM_ITERATIONS=10

mkdir -p /lizardfs/guarracino/chromosome_communities/boostraps
cd /lizardfs/guarracino/chromosome_communities/boostraps

python3 /lizardfs/guarracino/chromosome_communities/scripts/sample_contigs.py \
  /lizardfs/guarracino/chromosome_communities/pq_contigs/chr13.vs.chm13.100kbps.pq_contigs.union.fa.gz.fai \
  /lizardfs/guarracino/chromosome_communities/pq_contigs/chr14.vs.chm13.100kbps.pq_contigs.union.fa.gz.fai \
  /lizardfs/guarracino/chromosome_communities/pq_contigs/chr15.vs.chm13.100kbps.pq_contigs.union.fa.gz.fai \
  /lizardfs/guarracino/chromosome_communities/pq_contigs/chr21.vs.chm13.100kbps.pq_contigs.union.fa.gz.fai \
  /lizardfs/guarracino/chromosome_communities/pq_contigs/chr22.vs.chm13.100kbps.pq_contigs.union.fa.gz.fai \
  ${NUM_ITERATIONS} | sort | uniq > bootstrap.${NUM_ITERATIONS}_sets.txt

# Check (with lots of sets there could be duplicates)
cat bootstrap.${NUM_ITERATIONS}_sets.txt | wc -l
```

Prepare VCF files:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/boostraps/vcfs

cat /lizardfs/guarracino/chromosome_communities/recombination_rate/ref_names.txt | while read REF_NAME; do
    PATH_VCF_SNPS=$OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.vcf.gz
    
    i=1
    cat bootstrap.${NUM_ITERATIONS}_sets.txt | while read contigs; do
      #echo $(cat $f | paste -s -d',')
          
      bcftools view --samples $contigs $PATH_VCF_SNPS |\
       bgzip -@ $THREADS > /lizardfs/guarracino/chromosome_communities/boostraps/vcfs/bootstrap.1000sets.$REF_NAME.$i.vcf.gz
      tabix $PATH_VCF_SNPS
      
      i=$((i+1))
    done  
done
```

In parallel:

```shell
cat /lizardfs/guarracino/chromosome_communities/recombination_rate/ref_names.txt | while read REF_NAME; do
    echo $REF_NAME
    
    PATH_REF_FASTA=/lizardfs/guarracino/chromosome_communities/assemblies/$(echo $REF_NAME | tr '#' '.').fa.gz

    LENGTH_OF_SEQ=$(cut -f 2 ${PATH_REF_FASTA}.fai)
    START_OF_SEQ=0
    END_OF_SEQ=$LENGTH_OF_SEQ
    
    NUM_OF_SEGS=$(echo $LENGTH_OF_SEQ / $SEG_LENGTH | bc)
    NUM_OF_SEGS_PLUS_1=$(echo $LENGTH_OF_SEQ / $SEG_LENGTH + 1 | bc)
    
    
    seq 1 1000 | while read contigs; do
        #echo $(cat $f | paste -s -d',')
              
        PATH_VCF_SNPS=/lizardfs/guarracino/chromosome_communities/boostraps/vcfs/bootstrap.1000sets.$REF_NAME.$i.vcf.gz
        
        PATH_VCF_SNPS=$OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.vcf.gz
        
        N=$(zgrep '^#CHROM' $PATH_VCF_SNPS -m 1 | cut -f 10- | tr '\t' '\n' | grep grch38 -v | wc -l)
        
        #Unused   
        #VCF_SNPS_NAME=$(basename $PATH_VCF_SNPS .vcf.gz)
        #zgrep '^#' $PATH_VCF_SNPS -v | cut -f 2 | sort -k 1n > $VCF_SNPS_NAME.all.positions.txt
       
        cd /scratch/
        TEMP_DIR=${REF_NAME}_s${SEG_LENGTH}-temp
        mkdir -p ${REF_NAME}_s${SEG_LENGTH}-recomb/${TEMP_DIR}
        cd ${REF_NAME}_s${SEG_LENGTH}-recomb
    
        seq 1 $NUM_OF_SEGS_PLUS_1 | parallel -j $THREADS bash $pathLDhatChunk {} \
          $SEG_LENGTH \
          $NUM_OF_SEGS_PLUS_1 \
          $END_OF_SEQ \
          $N \
          $PATH_REF_FASTA \
          $REF_NAME \
          $TEMP_DIR \
          $PATH_VCF_SNPS \
          $runSumsLDhat \
          $pathLDhat \
          $runPhi \
          $pathGetIndexesR
        done
    done
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

  sample_with_genotypes=$(bcftools stats -s '-' $VCF_NAME.$chr.snps.tmp.vcf.gz | awk '$1=="PSC" && $12+$13>0 {print $3}' | paste -s -d',')
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
