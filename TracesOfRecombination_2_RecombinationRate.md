## Recombination rate

```shell
#pathLDhat=~/tools/LDhat
#runPhi=~/tools/PhiPack/src/Phi
#runSumsLDhat=/lizardfs/guarracino/chromosome_communities/scripts/Sums_LDHat_pack.sh
#pathGetIndexesR=/lizardfs/guarracino/chromosome_communities/scripts/get_indexes.R
#pathGetRecombinationRatePlotR=/lizardfs/guarracino/chromosome_communities/scripts/get_recombination_rate_plot.R

pathLDhat=~/git/LDhat/
runPhi=~/Downloads/Pangenomics/LDhat/PhiPack/src/Phi
runSumsLDhat=~/git/chromosome_communities/scripts/Sums_LDHat_pack.sh
pathGetIndexesR=~/git/chromosome_communities/scripts/get_indexes.R
pathGetRecombinationRatePlotR=~/git/chromosome_communities/scripts/get_recombination_rate_plot.R
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
#OUTPUT_DIR=/lizardfs/guarracino/chromosome_communities/recombination
OUTPUT_DIR=~/Downloads/Pangenomics/LDhat/recombination_rate

#PATH_VCF=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.100kbps.pq_contigs.union.s100k.l300k.p98.n97/chrACRO+refs.100kbps.pq_contigs.union.fa.gz.20c4357.4030258.41cabb1.smooth.fix.chm13.vcf.gz
PATH_VCF=~/Downloads/Pangenomics/LDhat/SimulatedPopulations/HatLandscapeN16Len1000000Nrhs15_th0.01_540_1.pggb/HatLandscapeN16Len1000000Nrhs15_th0.01_540_1.fixed.fa.5c3c9a3.7bdde5a.a933754.smooth.fix.gfa.vcf.gz
VCF_NAME=$(basename $PATH_VCF .vcf.gz)


# guix install vcftools
zgrep '^#' $PATH_VCF -v | cut -f 1 | sort | uniq | while read REF_NAME; do
  echo $REF_NAME
  
  #PATH_REF_FASTA=/lizardfs/guarracino/chromosome_communities/assemblies/chm13.$REF_NAME.fa.gz
  PATH_REF_FASTA=~/Downloads/Pangenomics/LDhat/SimulatedPopulations/HatLandscapeN16Len1000000Nrhs15_th0.01_540_1.fixed.$REF_NAME.fa

  PATH_VCF_SNPS=$OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.vcf.gz

  # Take SNPs and normalize
  vcftools --gzvcf $PATH_VCF --chr ${REF_NAME} --remove-indels --recode --recode-INFO-all --stdout |\
    bcftools norm -f $PATH_REF_FASTA -c s -m - |\
    bgzip -@ 48 > $OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.tmp.vcf.gz
  rm out.log
  
  # Take only samples with at least 1 variants in the VCF file
  sample_with_variants=$(bcftools stats -s '-' $OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.tmp.vcf.gz | awk '$1=="PSC" && $12+$13>0 {print $3}' | paste -s -d',')
  bcftools view --samples $sample_with_variants $OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.tmp.vcf.gz |\
   bgzip -@ 48 > $PATH_VCF_SNPS
  tabix $PATH_VCF_SNPS

  rm $OUTPUT_DIR/$VCF_NAME.$REF_NAME.snps.tmp.vcf.gz
done
```

Prepare ranges:

```shell
SEG_LENGTH=1000 # 1kb is recommended for LDJump

zgrep '^#' $PATH_VCF -v | cut -f 1 | sort | uniq | while read REF_NAME; do
    echo $REF_NAME
    
    #PATH_REF_FASTA=/lizardfs/guarracino/chromosome_communities/assemblies/chm13.$REF_NAME.fa.gz
    PATH_REF_FASTA=~/Downloads/Pangenomics/LDhat/SimulatedPopulations/HatLandscapeN16Len1000000Nrhs15_th0.01_540_1.fixed.$REF_NAME.fa
    
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
    
    
    TEMP_DIR=$OUTPUT_DIR/$REF_NAME-temp
    INDEXES_DIR=$OUTPUT_DIR/$REF_NAME-indexes
    
    mkdir -p $TEMP_DIR
    mkdir -p $INDEXES_DIR
    
    y=1
    seq 1 $NUM_OF_SEGS_PLUS_1 | while read x; do
      # Compute range
      ix=$(echo "($x - 1) * $SEG_LENGTH + $y" | bc)
      ex=$(echo "$ix + $SEG_LENGTH - $y" | bc)
      y=1
      fa_start=$(echo "($x - 1) * $SEG_LENGTH + $y" | bc)
      fa_end=$(echo "$fa_start + $SEG_LENGTH - $y" | bc)
      if [ $x == $NUM_OF_SEGS_PLUS_1 ]; then
        ix=$(echo "($x - 1) * $SEG_LENGTH + $y" | bc)
        ex=$endofseq
        fa_start=$(echo "($x - 1) * $SEG_LENGTH + $y" | bc)
        fa_end=$(echo "$fa_start + ($LENGTH_OF_SEQ - $fa_start)" | bc)
      fi;
      
      # Check if the range makes sense
      xx=$(echo "$ex - $ix" | bc)
      if (( $xx < 1 )); then
        continue
      fi
      
      # Check if there could be variants in the range
      if (( $ex < $MIN_POS )); then
        continue
      fi
        if (( $ix > $MAX_POS )); then
        continue
      fi
      
      
      # Extract reference in the range
      echo $x $ix $ex
      samtools faidx $PATH_REF_FASTA ${REF_NAME}:$fa_start-$fa_end > ${TEMP_DIR}/${REF_NAME}_$ix-$ex.fa
    
      # Extract variants in the range
      vcftools --gzvcf $PATH_VCF_SNPS --chr $REF_NAME \
        --from-bp $ix --to-bp $ex \
        --recode --recode-INFO-all \
        --stdout | bgzip -@ 16 -c > ${TEMP_DIR}/sel_${ix}_${ex}.vcf.gz
      tabix ${TEMP_DIR}/sel_${ix}_${ex}.vcf.gz
    
      num_variants=$(zgrep '^#' ${TEMP_DIR}/sel_${ix}_${ex}.vcf.gz -vc)
    
      if [ $num_variants == 0 ]; then
        echo -e "Segregating sites=0\nAverage PWD=0\nWatterson theta=0\nTajima D statistic=0\nFu and Li D* statistic=0\nVariance PWD=0\n" > Sums_part_main_job$x.txt 
      else
          # Generate the FASTA chunk
          zgrep '^#CHROM' ${TEMP_DIR}/sel_${ix}_${ex}.vcf.gz -m 1 | cut -f 10- | tr '\t' '\n' | while read SAMPLE; do
            echo $SAMPLE
        
            bcftools consensus -s $SAMPLE \
              -f ${TEMP_DIR}/${REF_NAME}_$ix-$ex.fa \
              ${TEMP_DIR}/sel_${ix}_${ex}.vcf.gz |\
              sed "s/${REF_NAME}/$SAMPLE/g" >> ${TEMP_DIR}/sel_${ix}_${ex}.fa
            #samtools faidx ${TEMP_DIR}/sel_${ix}_${ex}.fa
          done
          
          lpart=$(echo "$ex - $ix + 1" | bc)
          # Output: Sums_part_main_job$x.txt
          $runSumsLDhat \
            ${TEMP_DIR}/sel_${ix}_${ex}.fa \
            $N \
            $lpart \
            $pathLDhat \
            $x \
            "job" \
            ${TEMP_DIR}
      fi
      
      # guix install r-ape
      # guix install r-pegas
      # guix install r-adegenet
      Rscript $pathGetIndexesR $x ${TEMP_DIR}/sel_${ix}_${ex}.fa $runPhi "job" ${TEMP_DIR} | tail -n 1 \
        > ${TEMP_DIR}/${REF_NAME}.indexes.$x.tsv
      
#      if (( $x > 100 )); then
#        break
#      fi
    done
    
    cat ${TEMP_DIR}/${REF_NAME}.indexes.*.tsv > ${INDEXES_DIR}/indexes.tsv
    mv Sums_part_main_job*.txt ${INDEXES_DIR}
    rm Phi*
    rm out.log
    
    Rscript $pathGetRecombinationRatePlotR ${INDEXES_DIR} $SEG_LENGTH $N $LENGTH_OF_SEQ "job" ${INDEXES_DIR}
done


```


```shell


path_indexes_tsv<-'temp.indexes.all.tsv'
segLength<-1000
nn<-15
ll<-1000000
out<-'job'




# guix install r-data-table
```


```r

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

