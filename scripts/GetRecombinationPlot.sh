#!/bin/bash

i=$1
SEG_LENGTH=$2
REF_NAME=$3
BASE_DIR=$4
pathGetRecombinationRatePlotR=$5

#echo "i $i"
#echo "SEG_LENGTH $SEG_LENGTH"
#echo "REF_NAME $REF_NAME"
#echo "BASE_DIR $BASE_DIR"
#echo "pathGetRecombinationRatePlotR $pathGetRecombinationRatePlotR"

echo "$REF_NAME: iteration $i"

INDEXES_DIR=$BASE_DIR/${REF_NAME}_s${SEG_LENGTH}_i$i-indexes

# guix install r-data-table
Rscript $pathGetRecombinationRatePlotR \
  ${INDEXES_DIR}/Indexes.tsv \
  ${INDEXES_DIR}/Sums_part_main_job.txt \
  ${INDEXES_DIR}/Information.txt \
  ${INDEXES_DIR}/${REF_NAME}_s${SEG_LENGTH}_i$i \
  "$REF_NAME - iteration $i"
