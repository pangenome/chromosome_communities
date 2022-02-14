#!/bin/bash

x=$1
SEG_LENGTH=$2
NUM_OF_SEGS_PLUS_1=$3
END_OF_SEQ=$4
N=$5
ALL_VCF_POSITIONS=$6
PATH_REF_FASTA=$7
REF_NAME=$8
TEMP_DIR=$9
PATH_VCF_SNPS=${10}
runSumsLDhat=${11}
pathLDhat=${12}
runPhi=${13}
pathGetIndexesR=${14}


#echo "x - $x"
#echo "SEG_LENGTH - $SEG_LENGTH"
#echo "NUM_OF_SEGS_PLUS_1 - $NUM_OF_SEGS_PLUS_1"
#echo "END_OF_SEQ - $END_OF_SEQ"
#echo "N - $N"
#echo "ALL_VCF_POSITIONS - $ALL_VCF_POSITIONS"
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

# Padded counter to be able to retrieve the files already in order (with ls *.txt, cat *.txt, etc...)
num_of_digits=$(echo "${#NUM_OF_SEGS_PLUS_1} + 1" | bc)
printf -v padded_x "%0${num_of_digits}d" $x

# Extract reference in the range
echo $x $ix $ex
samtools faidx "$PATH_REF_FASTA" "${REF_NAME}":$fa_start-$fa_end > "${TEMP_DIR}"/"${REF_NAME}"_$ix-$ex.${padded_x}.fa

# Extract variants in the range
vcftools --gzvcf "$PATH_VCF_SNPS" --chr "$REF_NAME" \
  --from-bp $ix --to-bp $ex \
  --recode --recode-INFO-all \
  --stdout | bgzip -@ 1 -c > "${TEMP_DIR}"/sel_${ix}_${ex}.${padded_x}.vcf.gz
tabix "${TEMP_DIR}"/sel_${ix}_${ex}.${padded_x}.vcf.gz

#It is possible to have variants, but no samples with that variants (all genotypes '.' or '0')
# so it is better to check the number of samples with genotype equals to '1'.
###num_variants=$(zgrep '^#' "${TEMP_DIR}"/sel_${ix}_${ex}.${padded_x}.vcf.gz -vc)
sample_with_alternate_variants=$(bcftools stats -s '-' "${TEMP_DIR}"/sel_${ix}_${ex}.${padded_x}.vcf.gz | awk '$1=="PSC" && $13>0 {print $3}' | wc -l)

#Output (TAB-separated) : Segregating sites, Average PWD, Watterson theta, Tajima D statistic, Fu and Li D* statistic, Variance PWD
if [ $sample_with_alternate_variants == 0 ]; then
  echo -e "0\t0\t0\t0\t0\t0\t" > "${TEMP_DIR}"/"${REF_NAME}".sums_part_main.${padded_x}.txt
  echo -e "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t" > "${TEMP_DIR}"/"${REF_NAME}".indexes.${padded_x}.tsv
else
    # Generate the FASTA chunk
    zgrep '^#CHROM' "${TEMP_DIR}"/sel_${ix}_${ex}.${padded_x}.vcf.gz -m 1 | cut -f 10- | tr '\t' '\n' | grep grch38 -v | while read SAMPLE; do
      #echo $SAMPLE

      bcftools consensus -s "$SAMPLE" \
        -f "${TEMP_DIR}"/"${REF_NAME}"_$ix-$ex.${padded_x}.fa \
        "${TEMP_DIR}"/sel_${ix}_${ex}.${padded_x}.vcf.gz |\
        sed "s/${REF_NAME}/$SAMPLE/g" >> "${TEMP_DIR}"/sel_${ix}_${ex}.${padded_x}.fa
      #samtools faidx ${TEMP_DIR}/sel_${ix}_${ex}.${padded_x}.fa
    done
    # If it is needed to add also the reference itself
    #cat ${TEMP_DIR}/${REF_NAME}_$ix-$ex.${padded_x}.fa >> ${TEMP_DIR}/sel_${ix}_${ex}.${padded_x}.fa

    lpart=$(echo "$ex - $ix + 1" | bc)
    # Output: ${TEMP_DIR}/${REF_NAME}.sums_part_main.$x.txt
    $runSumsLDhat \
      "${TEMP_DIR}"/sel_${ix}_${ex}.${padded_x}.fa \
      $N \
      $lpart \
      "$pathLDhat" \
      ${padded_x} \
      "${REF_NAME}" \
      "${TEMP_DIR}"

    # guix install r-ape
    # guix install r-pegas
    # guix install r-adegenet
    Rscript "$pathGetIndexesR" ${padded_x} "${TEMP_DIR}"/sel_${ix}_${ex}.${padded_x}.fa "$runPhi" "${REF_NAME}" "${TEMP_DIR}" | tail -n 1 \
      > "${TEMP_DIR}"/"${REF_NAME}".indexes.${padded_x}.tsv
fi

#cat $ALL_VCF_POSITIONS | while read pos; do
#  if (($pos >= $ix && $pos <= $ex)); then
#
#
#    break
#  fi
#done
