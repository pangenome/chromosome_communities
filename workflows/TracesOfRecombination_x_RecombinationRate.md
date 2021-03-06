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

Choose the VCF file:

```shell
PATH_VCF=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.100kbps.pq_contigs.union.s100k.l300k.p98.n97/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.chm13.vcf.gz
#PATH_VCF=~/Downloads/Pangenomics/LDhat/SimulatedPopulations/HatLandscapeN16Len1000000Nrhs15_th0.01_540_1.pggb/HatLandscapeN16Len1000000Nrhs15_th0.01_540_1.fixed.fa.5c3c9a3.7bdde5a.a933754.smooth.fix.gfa.vcf.gz
#PATH_VCF=/lizardfs/guarracino/chromosome_communities/HatLandscapeN16Len1000000Nrhs15_th0.01_540_1.fixed.fa.5c3c9a3.7bdde5a.a933754.smooth.fix.gfa.vcf.gz

VCF_NAME=$(basename $PATH_VCF .vcf.gz)
SEG_LENGTH=10000 # 1kb is recommended for LDJump

# Get reference names
zgrep '^#' $PATH_VCF -v | cut -f 1 | sort | uniq > ref_names.txt
```

Prepare the VCF file with only SNPs:

```shell
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

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
  sample_with_genotypes=$(bcftools stats -s '-' $OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.tmp.vcf.gz | awk '$1=="PSC" && $12+$13>0 {print $3}' | paste -s -d',')
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

mkdir -p /lizardfs/guarracino/chromosome_communities/bootstraps
cd /lizardfs/guarracino/chromosome_communities/bootstraps

# Take contigs that are in the VCF files (that is, that have at least one valid genotype, 0 or 1)
cat /lizardfs/guarracino/chromosome_communities/recombination_rate/ref_names.txt | while read REF_NAME; do
  echo $REF_NAME

  PATH_VCF_SNPS=$OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.vcf.gz
  
  cat /lizardfs/guarracino/chromosome_communities/recombination_rate/ref_names.txt | while read REF_NAME_CONTIGS; do
    comm -12 \
      <(zgrep '^#CHROM' $PATH_VCF_SNPS -m 1 | cut -f 10- | tr '\t' '\n' | grep grch38 -v) \
      <(cut -f 1 /lizardfs/guarracino/chromosome_communities/pq_contigs/$(echo $REF_NAME_CONTIGS | cut -f 2 -d '#').vs.chm13.100kbps.pq_contigs.union.fa.gz.fai) \
      > $VCF_NAME.$REF_NAME.snp.contigs.$REF_NAME_CONTIGS.txt
  done
done

# Sample contigs
cat /lizardfs/guarracino/chromosome_communities/recombination_rate/ref_names.txt | while read REF_NAME; do
  echo $REF_NAME

  python3 /lizardfs/guarracino/chromosome_communities/scripts/sample_contigs.py \
    chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.chm13.$REF_NAME.snp.contigs.chm13#chr13.txt \
    chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.chm13.$REF_NAME.snp.contigs.chm13#chr14.txt \
    chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.chm13.$REF_NAME.snp.contigs.chm13#chr15.txt \
    chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.chm13.$REF_NAME.snp.contigs.chm13#chr21.txt \
    chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.chm13.$REF_NAME.snp.contigs.chm13#chr22.txt \
    ${NUM_ITERATIONS} | sort | uniq > bootstrap.$REF_NAME.${NUM_ITERATIONS}_sets.txt
done

# Check (with lots of sets there could be duplicates)
ls bootstrap*.txt | while read f; do num=$(cat $f | wc -l); echo "$f $num"; done
```

Prepare VCF files:

```shell
mkdir -p /lizardfs/guarracino/chromosome_communities/bootstraps/vcfs

cat /lizardfs/guarracino/chromosome_communities/recombination_rate/ref_names.txt | while read REF_NAME; do
    PATH_VCF_SNPS=$OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.vcf.gz
    
    i=1
    cat bootstrap.$REF_NAME.${NUM_ITERATIONS}_sets.txt | while read contigs; do
      echo "$REF_NAME - $i --- $contigs"
      
      PATH_VCF_SNPS_REF_SET_I=/lizardfs/guarracino/chromosome_communities/bootstraps/vcfs/bootstrap.$REF_NAME.${NUM_ITERATIONS}_sets.$i.vcf.gz
      bcftools view --samples $contigs $PATH_VCF_SNPS |\
       bgzip -@ $THREADS > $PATH_VCF_SNPS_REF_SET_I
      tabix $PATH_VCF_SNPS_REF_SET_I
      
      i=$((i+1))
    done  
done
```

In parallel:

```shell
cat /lizardfs/guarracino/chromosome_communities/recombination_rate/ref_names.txt | while read REF_NAME; do
    #echo $REF_NAME
    
    PATH_REF_FASTA=/lizardfs/guarracino/chromosome_communities/assemblies/$(echo $REF_NAME | tr '#' '.').fa.gz

    LENGTH_OF_SEQ=$(cut -f 2 ${PATH_REF_FASTA}.fai)
    START_OF_SEQ=0
    END_OF_SEQ=$LENGTH_OF_SEQ
    
    NUM_OF_SEGS=$(echo $LENGTH_OF_SEQ / $SEG_LENGTH | bc)
    NUM_OF_SEGS_PLUS_1=$(echo $LENGTH_OF_SEQ / $SEG_LENGTH + 1 | bc)
    
    seq 1 ${NUM_ITERATIONS} | while read i; do
        echo "$REF_NAME: iteration $i"
        PATH_VCF_SNPS_REF_SET_I=/lizardfs/guarracino/chromosome_communities/bootstraps/vcfs/bootstrap.$REF_NAME.${NUM_ITERATIONS}_sets.$i.vcf.gz
        
        #N=$(zgrep '^#CHROM' $PATH_VCF_SNPS_REF_SET_I -m 1 | cut -f 10- | tr '\t' '\n' | grep grch38 -v | wc -l)
        N=5 # We are sampling 5 contigs for each set

        #Unused   
        #VCF_SNPS_NAME=$(basename $PATH_VCF_SNPS .vcf.gz)
        #zgrep '^#' $PATH_VCF_SNPS -v | cut -f 2 | sort -k 1n > $VCF_SNPS_NAME.all.positions.txt
       
        cd /scratch/
        TEMP_DIR=${REF_NAME}_s${SEG_LENGTH}_i$i-temp
        mkdir -p ${REF_NAME}_s${SEG_LENGTH}_i$i-recomb/${TEMP_DIR}
        cd ${REF_NAME}_s${SEG_LENGTH}_i$i-recomb
    
        seq 1 $NUM_OF_SEGS_PLUS_1 | parallel -j $THREADS bash $pathLDhatChunk {} \
          $SEG_LENGTH \
          $NUM_OF_SEGS_PLUS_1 \
          $END_OF_SEQ \
          $N \
          $PATH_REF_FASTA \
          $REF_NAME \
          $TEMP_DIR \
          $PATH_VCF_SNPS_REF_SET_I \
          $runSumsLDhat \
          $pathLDhat \
          $runPhi \
          $pathGetIndexesR
        done
    done
done


```

```shell
BOOTSTRAP_OUTPUT_DIR=/lizardfs/guarracino/chromosome_communities/bootstraps

cat /lizardfs/guarracino/chromosome_communities/recombination_rate/ref_names.txt | while read REF_NAME; do
    #echo $REF_NAME
    
    seq 1 ${NUM_ITERATIONS} | while read i; do
        echo "$REF_NAME: iteration $i"

        cd /scratch/
        TEMP_DIR=${REF_NAME}_s${SEG_LENGTH}_i$i-recomb/${REF_NAME}_s${SEG_LENGTH}_i$i-temp
        
        INDEXES_DIR=${BOOTSTRAP_OUTPUT_DIR}/${REF_NAME}_s${SEG_LENGTH}_i$i-indexes
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
        N=5 # We are sampling 5 contigs for each set
        LENGTH_OF_SEQ=$(cut -f 2 ${PATH_REF_FASTA}.fai)
        printf $SEG_LENGTH"\t"$N"\t"$LENGTH_OF_SEQ"\n" > ${INDEXES_DIR}/Information.txt
    done
done
```


Recombination plots:

```shell
#scp -r guarracino@octopus04:/lizardfs/guarracino/chromosome_communities/bootstraps/chm13#chr*_s*_i*-indexes /home/guarracino/Downloads/Pangenomics/LDhat/bootstraps/

pathGetRecombinationRatePlotSh=/home/guarracino/git/chromosome_communities/scripts/GetRecombinationPlot.sh

#BASE_DIR=$BOOTSTRAP_OUTPUT_DIR
BASE_DIR=/home/guarracino/Downloads/Pangenomics/LDhat/bootstraps

cat ../recombination_rate/ref_names.txt | while read REF_NAME; do
    seq 1 ${NUM_ITERATIONS} | parallel -j $THREADS bash $pathGetRecombinationRatePlotSh {} \
      $SEG_LENGTH \
      $REF_NAME \
      $BASE_DIR \
      $pathGetRecombinationRatePlotR
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
  bcftools view --samples $sample_with_genotypes $VCF_NAME.$chr.snps.tmp.vcf.gz | bgzip -@ 48 > $VCF_NAME.$chr.snps.vcf.gz
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



[comment]: <> (```shell)

[comment]: <> (perl ~/git/vcf-conversion-tools/vcf2MS.pl example.vcf example.ms 629)

[comment]: <> (perl ~/git/vcf-conversion-tools/MS2LDhat.pl example.ms example.ldhat 629)

[comment]: <> (#or)

[comment]: <> (perl ~/git/vcf-conversion-tools/vcf2MS.pl example.10samples.vcf example.ms 11)

[comment]: <> (perl ~/git/vcf-conversion-tools/MS2LDhat.pl example.ms example 11)

[comment]: <> (~/git/LDhat/pairwise -seq example.ldhat.sites -loc example.ldhat.locs)

[comment]: <> (pip3 install biopython)

[comment]: <> (pip3 install pyfaidx)

[comment]: <> (pip3 install PyVCF)

[comment]: <> (python3 vcf2alignedFasta.py example.10samples.vcf.gz chr17.fasta )

[comment]: <> (# hack output &#40;num seq, length seq, 1 &#40;haploid&#41; / 2 &#40;diploid&#41;)

[comment]: <> (~/git/LDhat/convert -seq example.fasta )

[comment]: <> (~/git/LDhat/pairwise -seq sites.txt -loc locs.txt )

[comment]: <> (path_vcf=example.10samples.vcf)

[comment]: <> (path_vcf=pggb.wgg.88.chm13.1-22+X.norm.max50.vcf)

[comment]: <> (vk phylo fasta ${path_vcf} > ${path_vcf}.tmp.fasta # todo to improve)

[comment]: <> (seq_num=$&#40;grep '^>' ${path_vcf}.tmp.fasta -c&#41;)

[comment]: <> (seq_length=$&#40;head ${path_vcf}.tmp.fasta -n 2 | tail -n 1 | awk '{print length&#40;$0&#41;}'&#41;)

[comment]: <> (# 1 &#40;haploid&#41; / 2 &#40;diploid&#41;)

[comment]: <> (cat <&#40;echo -e $seq_num"\t"$seq_length"\t"1&#41; <&#40;cat ${path_vcf}.tmp.fasta&#41; > ${path_vcf}.fasta)

[comment]: <> (rm ${path_vcf}.tmp.fasta)

[comment]: <> (~/git/LDhat/convert -seq ${path_vcf}.fasta)

[comment]: <> (~/git/LDhat/pairwise -seq sites.txt -loc locs.txt )

[comment]: <> (```)
