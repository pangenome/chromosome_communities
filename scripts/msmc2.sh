#!/bin/sh

PATH_POPULATION_TXT=$1
PATH_HACKED_GENERATE_MULTIHETSEP_PY=$3
PATH_MSMC2=$3
PREFIX=$4

nsamples=$(cat $PATH_POPULATION_TXT | wc -l)
nsamples_minus_1=$(echo "$nsamples - 1" | bc)

samples=$(cat $PATH_POPULATION_TXT | tr '\n' ' ' | sed 's/.$//')

# Generate files for MSMC
# Using a hacked version, changing a few rows in generate_multihetsep.py to get phased fake-diploid data
#print ("HACKED: Non-diploid SNP found and considered as phased data: %s" % geno, file=sys.stderr)
#phased = True
#geno = "%s|%s" % (geno[0], geno[0])

# '#' breaks MSMC
python3 $PATH_HACKED_GENERATE_MULTIHETSEP_PY $samples | sed 's/#/-/g' > $PREFIX.multihetsep

# Estimating the effective population size
haplos=`seq 0 $nsamples_minus_1 | tr '\n' ',' | sed 's/.$//'`
\time -v $PATH_MSMC2 \
  -p 1*2+15*1+1*2 \
  -o $PREFIX.pop1 \
  -I $haplos \
  $PREFIX.multihetsep
