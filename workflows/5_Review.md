# Repetitive region centered on the SST1 array (chr13, chr14, chr21)

Flip the acrocentric graph:

```shell
RUN_ODGI=/home/guarracino/tools/odgi/bin/odgi-0b21b3525fc2ed9305c7df2386475a008a9337bd 

path_input_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.og
prefix=$(basename $path_input_og .og)

sbatch -p workers -c 48 --job-name acro-flip --wrap "hostname; cd /scratch && $RUN_ODGI flip -i $path_input_og -o ${prefix}.flip.og -t 48 -P; mv ${prefix}.flip.og /lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/"
```

Save the path names with the information of the flipped path, and rename them by removing the '_inv' suffix:

```shell
path_flipped_og=/lizardfs/guarracino/chromosome_communities/graphs/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.s50k.l250k.p98.n162/chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.flip.og
prefix=$(basename $path_flipped_og .og)

$RUN_ODGI paths -i $path_flipped_og -L > $prefix.path_names.txt

$RUN_ODGI view -i $path_flipped_og -g | sed '/^P/ s/_inv//' > $prefix.gfa
rm $path_flipped_og
$RUN_ODGI build -g $prefix.gfa -o $path_flipped_og -t 48 -P
```

Extract the region:

```shell
grep censat /lizardfs/guarracino/chromosome_communities/data/annotation/chm13.SST1.200kbps.aproximate.acros.bed | \
  bedtools slop -i - -b 1000000 -g <(cut -f 1,3 /lizardfs/guarracino/chromosome_communities/data/annotation/chm13.q_arms.approximate.acros.bed) \
  > SST1_regions_1Mpbs.bed

path_input_og=chrACRO+refs.pq_contigs.1kbps.hg002prox.hg002hifi.fa.gz.7ef1ba2.04f1c29.ebc49e1.smooth.final.og
$RUN_ODGI extract \
  -i $path_input_og -b SST1_regions_1Mpbs.bed \
  -d 100000 -e 3 -o - -t 48 -P | \
  $RUN_ODGI sort -i - -o SST1_regions_1Mpbs.og -O -p gYs -t 48 -P

# All targets are flipped in the flipped graph, so we flip p/q-arms coordinates for filtering
rm SST1_regions_1Mpbs.flip.bed
for ref in chm13#chr13 chm13#chr14 chm13#chr21; do
  len=$(grep $ref /lizardfs/guarracino/chromosome_communities/pq_contigs/chrACRO+refs.pq_contigs.1kbps.hg002prox.fa.gz.fai | cut -f 2)
  grep censat /lizardfs/guarracino/chromosome_communities/data/annotation/chm13.SST1.200kbps.aproximate.acros.bed | grep $ref | \
    bedtools slop -i - -b 1000000 -g <(cut -f 1,3 /lizardfs/guarracino/chromosome_communities/data/annotation/chm13.q_arms.approximate.acros.bed) |
    awk -v OFS='\t' -v len=$len '{print($1,len-$3,len-$2)}' | bedtools sort \
    >> SST1_regions_1Mpbs.flip.bed
done
$RUN_ODGI extract \
  -i $path_flipped_og -b SST1_regions_1Mpbs.bed \
  -d 100000 -e 3 -o - -t 48 -P | \
  $RUN_ODGI sort -i - -o SST1_regions_1Mpbs.flip.og -O -p gYs -t 48 -P
  
$RUN_ODGI view -i SST1_regions_2Mpbs.og -g > SST1_regions_2Mpbs.gfa
$RUN_ODGI layout -i SST1_regions_2Mpbs.og -T SST1_regions_2Mpbs.layout.tsv -t 48 -P

~/git/gfaestus/target/release/gfaestus SST1_regions_2Mpbs.gfa SST1_regions_2Mpbs.layout.tsv 



$RUN_ODGI extract \
  -i $path_input_og -b /lizardfs/guarracino/chromosome_communities/data/annotation/chm13.p_arms.approximate.acros.bed \
  -d 100000 -e 3 -o - -t 48 -P | \
  $RUN_ODGI sort -i - -o p_arms.og -O -p gYs -t 48 -P
```
