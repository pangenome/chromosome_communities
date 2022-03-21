#!/bin/sh

PATH_POPULATION_1_TXT=$1
PATH_POPULATION_2_TXT=$2
PATH_HACKED_GENERATE_MULTIHETSEP_PY=$3
PATH_COMBINE_CROSS_COAL_PY=$4
PATH_MSMC2=$5
PREFIX=$6

# It assumes population of the same size
nsamples=$(cat $PATH_POPULATION_1_TXT | wc -l)
nsamples_minus_1=$(echo "$nsamples - 1" | bc)
nsamples_mul_2_minus_1=$(echo "$nsamples * 2 - 1" | bc)

samples=$(cat $PATH_POPULATION_1_TXT $PATH_POPULATION_2_TXT | tr '\n' ' ' | sed 's/.$//'`)

# Generate files for MSMC
# Using a hacked version, changing a few rows in generate_multihetsep.py to get phased fake-diploid data
#print ("HACKED: Non-diploid SNP found and considered as phased data: %s" % geno, file=sys.stderr)
#phased = True
#geno = "%s|%s" % (geno[0], geno[0])
python3 $PATH_HACKED_GENERATE_MULTIHETSEP_PY $samples > $PREFIX.multihetsep

# Estimating the effective population size
haplos1=`seq 0 $nsamples_minus_1 | tr '\n' ',' | sed 's/.$//'`
\time -v $PATH_MSMC2 \
  -p 1*2+15*1+1*2 \
  -o $PREFIX.pop1 \
  -I $haplos1 \
  $PREFIX.multihetsep

haplos2=`seq $nsamples $nsamples_mul_2_minus_1 | tr '\n' ',' | sed 's/.$//'`
\time -v $PATH_MSMC2 \
  -p 1*2+15*1+1*2 \
  -o $PREFIX.pop2 \
  -I $haplos2 \
  $PREFIX.multihetsep


# Estimating population separation history
combinations=$((for x in `seq 0 $nsamples_minus_1`; do
  for y in `seq $nsamples $nsamples_mul_2_minus_1`; do
    echo "$x-$y"
  done
done) | tr '\n' ',' | sed 's/.$//')

# --skipAmbiguous: skip sites with ambiguous phasing. Recommended for cross population analysis
\time -v $PATH_MSMC2 \
  -p 1*2+15*1+1*2 \
  --skipAmbiguous \
  -o $PREFIX.pop1and2 \
  -I $combinations \
  $PREFIX.multihetsep


# Since we have obtained the coalescence rates independently, we have allowed MSMC2 to choose different
# time interval boundaries in each case, depending on the observed heterozygosity within and across populations.
# We therefore combine the three rates to obtain cross-population results.
python3 $PATH_COMBINE_CROSS_COAL_PY \
  $PREFIX.pop1and2.final.txt \
  $PREFIX.pop1.final.txt \
  $PREFIX.pop2.final.txt \
  > $PREFIX.pop1and2.combined.final.txt
