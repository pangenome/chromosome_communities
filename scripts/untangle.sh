path_input_og_gz=$1
references=$2
chromosome=$3
merge_dist=$4
prefix=$5

pattern_references=$(echo "$references" | sed 's/ /\\|/')
pattern_chromosomes=$(echo "$chromosome" | sed 's/ /\\|/')

suffix_references=$(echo "$references" | sed 's/ /_/')

\time -v odgi untangle \
  -i <(zcat $path_input_og_gz) \
  -R <(odgi paths -i <(zcat $path_input_og_gz) -L | grep $pattern_references | grep $pattern_chromosomes) \
  -m $merge_dist \
  -s 0 -j 0 \
  --paf-output \
  -t 48 -P > $prefix.all.vs.$suffix_references.m$merge_dist.paf
