#!/bin/bash

x=$1
SEG_LENGTH=$2
NUM_OF_SEGS_PLUS_1=$3
END_OF_SEQ=$4
N=$5
PATH_REF_FASTA=$6
REF_NAME=$7
TEMP_DIR=$8
PATH_VCF_SNPS=$9
runSumsLDhat=${10}
pathLDhat=${11}
runPhi=${12}
pathGetIndexesR=${13}


#echo "x - $x"
#echo "SEG_LENGTH - $SEG_LENGTH"
#echo "NUM_OF_SEGS_PLUS_1 - $NUM_OF_SEGS_PLUS_1"
#echo "END_OF_SEQ - $END_OF_SEQ"
#echo "N - $N"
#echo "PATH_REF_FASTA - $PATH_REF_FASTA"
#echo "REF_NAME - $REF_NAME"
#echo "TEMP_DIR - $TEMP_DIR"
#echo "PATH_VCF_SNPS - $PATH_VCF_SNPS"
#echo "runSumsLDhat - $runSumsLDhat"
#echo "pathLDhat - $pathLDhat"
#echo "runPhi - $runPhi"
#echo "pathGetIndexesR - $pathGetIndexesR"

y=1
# Compute range
ix=$(echo "($x - 1) * $SEG_LENGTH + $y" | bc)
ex=$(echo "$ix + $SEG_LENGTH - $y" | bc)
y=1
fa_start=$(echo "($x - 1) * $SEG_LENGTH + $y" | bc)
fa_end=$(echo "$fa_start + $SEG_LENGTH - $y" | bc)
if [ "$x" == "$NUM_OF_SEGS_PLUS_1" ]; then
  ix=$(echo "($x - 1) * $SEG_LENGTH + $y" | bc)
  ex=$END_OF_SEQ
  fa_start=$(echo "($x - 1) * $SEG_LENGTH + $y" | bc)
  fa_end=$(echo "$fa_start + ($LENGTH_OF_SEQ - $fa_start)" | bc)
fi;

# Check if the range makes sense
xx=$(echo "$ex - $ix" | bc)
if (( $xx < 1 )); then
  exit
fi

# Extract reference in the range
echo $x $ix $ex
samtools faidx $PATH_REF_FASTA ${REF_NAME}:$fa_start-$fa_end > ${TEMP_DIR}/${REF_NAME}_$ix-$ex.fa

# Extract variants in the range
vcftools --gzvcf $PATH_VCF_SNPS --chr $REF_NAME \
  --from-bp $ix --to-bp $ex \
  --recode --recode-INFO-all \
  --stdout | bgzip -@ 1 -c > ${TEMP_DIR}/sel_${ix}_${ex}.vcf.gz
tabix ${TEMP_DIR}/sel_${ix}_${ex}.vcf.gz

num_variants=$(zgrep '^#' ${TEMP_DIR}/sel_${ix}_${ex}.vcf.gz -vc)

# Padded counter to be able to retrieve the files already in order (with ls *.txt, cat *.txt, etc...)
num_of_digits=$(echo "${#NUM_OF_SEGS_PLUS_1} + 1" | bc)
printf -v padded_x "%0${num_of_digits}d" $x

#Output (TAB-separated) : Segregating sites, Average PWD, Watterson theta, Tajima D statistic, Fu and Li D* statistic, Variance PWD
if [ $num_variants == 0 ]; then
  echo -e "0\t0\t0\t0\t0\t0\t" > Sums_part_main_job${padded_x}.txt
else
    # Generate the FASTA chunk
    zgrep '^#CHROM' ${TEMP_DIR}/sel_${ix}_${ex}.vcf.gz -m 1 | cut -f 10- | tr '\t' '\n' | grep grch38 -v | while read SAMPLE; do
      #echo $SAMPLE

      bcftools consensus -s $SAMPLE \
        -f ${TEMP_DIR}/${REF_NAME}_$ix-$ex.fa \
        ${TEMP_DIR}/sel_${ix}_${ex}.vcf.gz |\
        sed "s/${REF_NAME}/$SAMPLE/g" >> ${TEMP_DIR}/sel_${ix}_${ex}.fa
      #samtools faidx ${TEMP_DIR}/sel_${ix}_${ex}.fa
    done
    # Add also the reference itself
    #cat ${TEMP_DIR}/${REF_NAME}_$ix-$ex.fa >> ${TEMP_DIR}/sel_${ix}_${ex}.fa

    lpart=$(echo "$ex - $ix + 1" | bc)
    # Output: Sums_part_main_job$x.txt
    $runSumsLDhat \
      ${TEMP_DIR}/sel_${ix}_${ex}.fa \
      $N \
      $lpart \
      $pathLDhat \
      ${padded_x} \
      "job" \
      ${TEMP_DIR}
fi

# guix install r-ape
# guix install r-pegas
# guix install r-adegenet
Rscript $pathGetIndexesR ${padded_x} ${TEMP_DIR}/sel_${ix}_${ex}.fa $runPhi "job" ${TEMP_DIR} | tail -n 1 \
  > ${TEMP_DIR}/${REF_NAME}.indexes.${padded_x}.tsv
